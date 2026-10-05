#= T21L20 moist BM validation with an independently observed physics budget.

Run from the repository root:
    julia --project=. exp/BettsMiller_Validation/validate.jl OUTPUT [DAYS=5] [DT=600]
Raw fields and hourly observations stay in OUTPUT; commit the small result
report, not generated model data. The production dynamics and physics kernels
advance the model. A separate replay observes individual physics processes and
must match the production physics state and diagnostics at every step.
=#
using JGCM, Printf, Serialization, SHA, TOML

const AP = JGCM.Atmos_Param_Module
const SD = JGCM.Spectral_Dynamics_Module
const PG = JGCM.Press_And_Geopot_Module
const TI = JGCM.Time_Integrator_Module
const SI = JGCM.Semi_Implicit_Module

prescribed_sst(longitude, latitude) = 271.0 + 29.0 * exp(-latitude^2 / (2 * (26pi / 180)^2))

function configuration(output, days, dt)
    duration = round(Int, days * 86400)
    params = Dict{String,Any}(
        "do_mass_correction" => true, "do_energy_correction" => true,
        "do_water_correction" => true, "use_virtual_temperature" => true,
        "do_Betts_Miller" => true, "bm_tau" => 7200.0,
        "bm_relative_humidity" => 0.8, "bm_energy_correction" => :isca,
        "do_Lscale_Cond" => true, "condensation_heating_fraction" => 1.0,
        "do_LRF" => false, "do_Sensible_Heating" => true, "C_H" => 0.0044,
        "do_Surface_Evaporation" => true, "C_E" => 0.0044,
        "lower_boundary_temperature" => prescribed_sst,
        "do_Implicit_PBL_Scheme" => true, "C_D" => 0.0044,
        "PBL_Top_Mode" => :PressureLevel, "PBL_Top_Value" => 85000.0,
        "do_HS_Forcing" => true, "σ_b" => 0.7, "k_a" => 1 / 40,
        "k_s" => 1 / 4, "k_f" => 1.0, "ΔT_y" => 60.0, "Δθ_z" => 10.0,
    )
    return Model_Config(
        name = "BM_validation_dt$dt", model_type = :PrimitiveEquation,
        num_fourier = 21, nθ = 32, nd = 20, radius = 6371.0e3,
        omega = 7.292e-5, grav = 9.8, vert_coord_option = "even_sigma",
        vert_difference_option = "simmons_and_burridge", vert_ref_level_option = "second_centered_wts",
        Δt = dt, end_time = duration, day_to_sec = 86400, damping_order = 4,
        damping_coef = 1.15741e-4, robert_coef = 0.04, implicit_coef = 0.5,
        initial_condition = :Moist_Spinup, moisture_processes = true,
        output_path = output, output_filename = joinpath(output, "unused.nc"),
        logger = joinpath(output, "progress.log"), vars_to_output = Symbol[],
        output_interval = 3600, saving_frequency = 0, physics_params = params,
    )
end

function frame(config)
    JGCM.Driver.validate_config(config)
    nlon, nlat, nd = 2config.nθ, config.nθ, config.nd
    mesh = Spectral_Spherical_Mesh(21, 22, nlon, nlat, nd, config.radius)
    vert = Vert_Coordinate(nlon, nlat, nd, config.vert_coord_option,
        config.vert_difference_option, config.vert_ref_level_option)
    atmo = Atmo_Data(config.name, nlon, nlat, nd, true, true, true, true, mesh.sinθ;
        radius = config.radius, omega = config.omega, grav = config.grav)
    integrator = TI.Filtered_Leapfrog(config.robert_coef, config.damping_order,
        config.damping_coef, mesh.laplacian_eig, config.implicit_coef, config.Δt,
        true, 0, config.end_time)
    semi = SI.Semi_Implicit_Solver(vert, atmo, integrator, 1.0e5, fill(300.0, nd), mesh.wave_numbers)
    dyn = Dyn_Data(config.name, 21, 22, nlon, nlat, nd)
    config.physics_params["BM_state"] = JGCM.Driver.betts_miller_state(config.physics_params, nd, config.Δt, true)
    Initialize_Atmos_State!(mesh, atmo, dyn, vert, config)
    return (; config, mesh, vert, atmo, integrator, semi, dyn)
end

"""Area mean [water kg/m², cp*T+Lv*q+K J/m², dry mass kg/m²]."""
function ledger(mesh, atmo, u, v, temperature, humidity, ph)
    water, energy, dry = 0.0, 0.0, 0.0
    @inbounds for k in axes(temperature, 3), j in axes(temperature, 2), i in axes(temperature, 1)
        mass_weight = (ph[i,j,k+1] - ph[i,j,k]) * mesh.wts[j] / (2mesh.nλ * atmo.grav)
        q = humidity[i,j,k]
        water += q * mass_weight
        energy += (atmo.cp_air * temperature[i,j,k] + atmo.Lv * q +
            0.5 * (u[i,j,k]^2 + v[i,j,k]^2)) * mass_weight
        dry += (1 - q) * mass_weight
    end
    return [water, energy, dry]
end

function interfaces(vert, ps)
    return reshape(vert.ak, 1, 1, :) .+ reshape(vert.bk, 1, 1, :) .* ps
end

function dyn_ledger(f, level)
    d = f.dyn
    u, v, t, q, ps = (getfield(d, Symbol("grid_", key, "_", level)) for key in (:u, :v, :t, :q, :ps))
    return ledger(f.mesh, f.atmo, u, v, t, q, interfaces(f.vert, ps))
end

function observer(f)
    shape = size(f.dyn.grid_t_c)
    nlon, nlat, nd = shape
    return (;
        work = AP.Physics_Workspace(nlon, nlat, nd), pf = similar(f.dyn.grid_p_full),
        ph = similar(f.dyn.grid_p_half), ps = similar(f.dyn.grid_ps_c),
        dp = similar(f.dyn.grid_Δp), logpf = similar(f.dyn.grid_p_full), logph = similar(f.dyn.grid_p_half),
        mixing = similar(f.dyn.K_E), teq = zeros(shape),
        average_bm_rain = zeros(nlon,nlat,1), average_rain = zeros(nlon,nlat,1),
        average_lh = zeros(nlon,nlat,1), average_sh = zeros(nlon,nlat,1),
        average_bm_t = zeros(shape), average_bm_q = zeros(shape),
        previous_regime = fill(-1, nlon, nlat),
    )
end

function replay!(f, o, statistics)
    config, mesh, vert, atmo, integrator, semi, dyn = f
    p, w = config.physics_params, o.work
    u, v, t, q = w.grid_u, w.grid_v, w.grid_t, w.grid_q
    for (target, source) in zip((u,v,t,q,o.ps,o.pf,o.ph,o.dp),
            (dyn.grid_u_n,dyn.grid_v_n,dyn.grid_t_n,dyn.grid_q_n,dyn.grid_ps_n,dyn.grid_p_full,dyn.grid_p_half,dyn.grid_Δp))
        copyto!(target, source)
    end
    for field in (o.average_bm_rain,o.average_rain,o.average_bm_t,o.average_bm_q,o.average_lh,o.average_sh)
        fill!(field, 0.0)
    end
    dt = config.Δt
    nsubsteps = TI.Get_Δt(integrator) ÷ dt
    stages = Dict{String,Vector{Float64}}()
    rates = Dict("bm_rain" => 0.0, "lscale_rain" => 0.0, "evaporation" => 0.0)
    before = ledger(mesh, atmo, u,v,t,q,o.ph)
    function observe(name)
        after = ledger(mesh, atmo, u,v,t,q,o.ph)
        increment = after - before
        stages[name] = get(stages, name, zeros(3)) + increment
        before = after
        return increment
    end
    for substep = 1:nsubsteps
        dry_at_substep_start = before[3]
        AP._reset_substep_diagnostics!(w)
        copyto!(w.grid_q_before, q)
        @. w.grid_Δp_before = o.ph[:,:,2:end] - o.ph[:,:,1:end-1]
        Betts_Miller!(p["BM_state"], atmo, t,q,o.pf,o.ph,
            w.grid_bm_t_tendency,w.grid_bm_q_tendency,w.grid_bm_precip)
        for j in axes(t,2), i in axes(t,1)
            regime = w.grid_bm_precip[i,j,1] > 0 ? 2 :
                any(!iszero, @view(w.grid_bm_t_tendency[i,j,:])) ||
                any(!iszero, @view(w.grid_bm_q_tendency[i,j,:])) ? 1 : 0
            statistics[("none_samples","shallow_samples","deep_samples")[regime+1]] += 1
            dtemperature = dt * maximum(abs, @view(w.grid_bm_t_tendency[i,j,:]))
            dhumidity = dt * maximum(abs, @view(w.grid_bm_q_tendency[i,j,:]))
            statistics["max_bm_temperature_increment"] = max(statistics["max_bm_temperature_increment"], dtemperature)
            statistics["max_bm_humidity_increment"] = max(statistics["max_bm_humidity_increment"], dhumidity)
            if o.previous_regime[i,j] >= 0 && o.previous_regime[i,j] != regime
                statistics["regime_changes"] += 1
                statistics["max_transition_temperature_increment"] = max(statistics["max_transition_temperature_increment"], dtemperature)
                statistics["max_transition_humidity_increment"] = max(statistics["max_transition_humidity_increment"], dhumidity)
            end
            o.previous_regime[i,j] = regime
        end
        @. t += dt * w.grid_bm_t_tendency
        @. q += dt * w.grid_bm_q_tendency
        copyto!(w.grid_precip, w.grid_bm_precip)
        delta = observe("betts_miller")
        bm_rain = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(mesh, w.grid_bm_precip)
        statistics["max_bm_water_residual"] = max(statistics["max_bm_water_residual"], abs(delta[1] + dt*bm_rain))
        statistics["max_bm_energy_residual"] = max(statistics["max_bm_energy_residual"], abs(delta[2]))
        abs(delta[1] + dt*bm_rain) < 1.0e-9 || error("BM water budget failed")
        abs(delta[2]) < 0.01 || error("BM fixed-mass energy budget failed")
        AP.Lscale_Cond!(atmo,t,q,o.pf,o.ph,dt,1.0,w.grid_lscale_t_tendency,
            w.grid_lscale_q_tendency,w.grid_liquid_water_content,w.grid_precip)
        delta = observe("large_scale_condensation")
        total_rain = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(mesh, w.grid_precip)
        statistics["max_lscale_water_residual"] = max(statistics["max_lscale_water_residual"], abs(delta[1] + dt*(total_rain-bm_rain)))
        statistics["max_lscale_energy_residual"] = max(statistics["max_lscale_energy_residual"], abs(delta[2]))
        abs(delta[1] + dt*(total_rain-bm_rain)) < 1.0e-9 || error("Condensation water budget failed")
        abs(delta[2]) < 0.01 || error("Condensation fixed-mass energy budget failed")
        wind, height = AP.Calculate_V_c_za!(w.pbl,atmo,o.ph,o.ps,u,v,t,q)
        AP.Sensible_Heating!(mesh,atmo,o.ph,t,w.grid_shflx,wind,height,dt,p["C_H"],prescribed_sst)
        delta = observe("sensible_heat")
        sensible = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(mesh,w.grid_shflx)
        statistics["max_sensible_energy_residual"] = max(statistics["max_sensible_energy_residual"],abs(delta[2]-dt*sensible))
        AP.Surface_Evaporation!(mesh,atmo,o.ps,o.ph,q,w.grid_lhflx,wind,height,dt,p["C_E"],prescribed_sst)
        delta = observe("surface_evaporation")
        evaporation = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(mesh,w.grid_lhflx) / atmo.Lv
        statistics["max_evaporation_water_residual"] = max(statistics["max_evaporation_water_residual"],abs(delta[1]-dt*evaporation))
        abs(delta[1]-dt*evaporation) < 1.0e-9 || error("Evaporation water budget failed")
        AP.Implicit_PBL_Mixing!(w.pbl,atmo,o.pf,o.ph,t,q,o.mixing,wind,height,p,dt,p["C_D"])
        observe("pbl_mixing")
        AP.Rayleigh_Friction!(atmo,dt,config.day_to_sec,o.ph,o.pf,u,v,t,p)
        observe("rayleigh_with_heating")
        AP.Newtonian_Relaxation!(atmo,dt,config.day_to_sec,mesh.sinθ,o.ph,o.pf,t,o.teq,p)
        observe("newtonian_relaxation")
        AP._validate_and_clean_physics_state!(u,v,t,q)
        observe("humidity_cleanup")
        AP.Dry_Air_Adjustment!(vert,o.ps,w.grid_q_before,q,w.grid_Δp_before)
        PG.Pressure_Variables!(vert,o.ps,o.ph,o.dp,o.logph,o.pf,o.logpf)
        delta = observe("pressure_adjustment")
        statistics["max_pressure_water_residual"] = max(statistics["max_pressure_water_residual"],abs(delta[1]))
        statistics["max_pressure_dry_mass_residual"] = max(statistics["max_pressure_dry_mass_residual"],abs(before[3] - dry_at_substep_start))
        abs(before[3] - dry_at_substep_start) < 1.0e-8 || error("Pressure adjustment dry-mass budget failed")
        for (average, rate) in zip((o.average_bm_rain,o.average_rain,o.average_bm_t,o.average_bm_q,o.average_lh,o.average_sh),
                (w.grid_bm_precip,w.grid_precip,w.grid_bm_t_tendency,w.grid_bm_q_tendency,w.grid_lhflx,w.grid_shflx))
            @. average += rate / nsubsteps
        end
        rates["bm_rain"] += bm_rain / nsubsteps
        rates["lscale_rain"] += (total_rain - bm_rain) / nsubsteps
        rates["evaporation"] += evaporation / nsubsteps
    end
    return stages, rates
end

function state_hash(dyn)
    io = IOBuffer()
    for key in (:u,:v,:t,:q,:ps)
        write(io, getfield(dyn,Symbol("grid_",key,"_c")))
    end
    return bytes2hex(sha256(take!(io)))
end

function run_validation(output, days, dt)
    mkpath(output)
    config = configuration(output,days,dt)
    f = frame(config)
    o = observer(f)
    stats = Dict{String,Float64}(key => 0.0 for key in (
        "none_samples","shallow_samples","deep_samples","regime_changes",
        "max_bm_temperature_increment","max_bm_humidity_increment",
        "max_transition_temperature_increment","max_transition_humidity_increment",
        "max_bm_water_residual","max_bm_energy_residual","max_lscale_water_residual",
        "max_lscale_energy_residual","max_sensible_energy_residual","max_evaporation_water_residual",
        "max_pressure_water_residual","max_pressure_dry_mass_residual","max_replay_error",
        "max_ledger_water_residual","max_ledger_energy_residual"))
    initial = dyn_ledger(f,:c)
    initial_hash = state_hash(f.dyn)
    serialize(joinpath(output,"initial.bin"),(u=copy(f.dyn.grid_u_c),v=copy(f.dyn.grid_v_c),
        t=copy(f.dyn.grid_t_c),q=copy(f.dyn.grid_q_c),ps=copy(f.dyn.grid_ps_c)))
    previous, current = Dict{String,Vector{Float64}}(), Dict{String,Vector{Float64}}()
    rain_previous, rain_current = zeros(3), zeros(3)
    bounds = [minimum(f.dyn.grid_t_c),maximum(f.dyn.grid_t_c),minimum(f.dyn.grid_q_c),maximum(f.dyn.grid_q_c)]
    hourly = open(joinpath(output,"hourly.tsv"),"w")
    println(hourly,"time_s\tmin_T\tmax_T\tmin_q\tmax_q\twater\tenergy\tbm_rain_mm\tlscale_rain_mm\tevaporation_mm\tbm_rate\tlscale_rate")
    steps = JGCM.Driver.time_steps(config.end_time,dt)
    start = time_ns()
    try
        for step = 1:steps
            base = dyn_ledger(f,:p)
            SD.Spectral_Dynamics!(config,f.mesh,f.vert,f.atmo,f.dyn,f.semi;advance_time=false)
            PG.Pressure_Variables!(f.vert,f.dyn.grid_ps_n,f.dyn.grid_p_half,f.dyn.grid_Δp,
                f.dyn.grid_lnp_half,f.dyn.grid_p_full,f.dyn.grid_lnp_full)
            dynamics_increment = dyn_ledger(f,:n) - base
            stages, rates = replay!(f,o,stats)
            stages["dynamics_with_corrections"] = dynamics_increment
            AP.Spectral_Physics!(config,f.mesh,f.vert,f.atmo,f.dyn,f.semi,config.physics_params)
            for (actual,expected) in zip((f.dyn.grid_u_n,f.dyn.grid_v_n,f.dyn.grid_t_n,f.dyn.grid_q_n,f.dyn.grid_ps_n,
                    f.dyn.grid_bm_t_tendency,f.dyn.grid_bm_q_tendency,f.dyn.grid_bm_precip,f.dyn.grid_precip,f.dyn.grid_lhflx,f.dyn.grid_shflx),
                    (o.work.grid_u,o.work.grid_v,o.work.grid_t,o.work.grid_q,o.ps,
                    o.average_bm_t,o.average_bm_q,o.average_bm_rain,o.average_rain,o.average_lh,o.average_sh))
                stats["max_replay_error"] = max(stats["max_replay_error"],maximum(abs.(actual-expected)))
                isapprox(actual,expected;rtol=2.0e-14,atol=1.0e-14) || error("Observed replay differs from production physics")
            end
            before = dyn_ledger(f,:n)
            SD._synchronize_physics_next!(f.mesh,f.vert,f.atmo,f.dyn)
            final = dyn_ledger(f,:n)
            stages["spectral_synchronization"] = final - before
            # Leapfrog n is built from p, not c. Keep the budget history in
            # the same recurrence instead of summing 2dt over every sample.
            next = Dict(name => get(previous,name,zeros(3)) + value for (name,value) in stages)
            predicted = initial + reduce(+,values(next))
            stats["max_ledger_water_residual"] = max(stats["max_ledger_water_residual"],abs(final[1]-predicted[1]))
            stats["max_ledger_energy_residual"] = max(stats["max_ledger_energy_residual"],abs(final[2]-predicted[2]))
            abs(final[1]-predicted[1]) < 1.0e-8 || error("Full-model water ledger failed")
            abs(final[2]-predicted[2]) < 0.1 || error("Full-model energy ledger failed")
            interval = TI.Get_Δt(f.integrator)
            rain_next = rain_previous + interval * [rates["bm_rain"],rates["lscale_rain"],rates["evaporation"]]
            previous, current = current, next
            rain_previous, rain_current = rain_current, rain_next
            JGCM.Dyn_Data_Module.Time_Advance!(f.dyn)
            PG.Compute_Pressures_And_Heights!(f.atmo,f.vert,f.dyn.grid_ps_c,f.dyn.grid_geopots,f.dyn.grid_t_c,
                f.dyn.grid_p_half,f.dyn.grid_Δp,f.dyn.grid_lnp_half,f.dyn.grid_p_full,f.dyn.grid_lnp_full,
                f.dyn.grid_z_full,f.dyn.grid_z_half,f.dyn.grid_q_c)
            f.integrator.init_step && SI.Update_Init_Step!(f.semi)
            f.integrator.time += dt
            extrema_now = [minimum(f.dyn.grid_t_c),maximum(f.dyn.grid_t_c),minimum(f.dyn.grid_q_c),maximum(f.dyn.grid_q_c)]
            bounds = [min(bounds[1],extrema_now[1]),max(bounds[2],extrema_now[2]),min(bounds[3],extrema_now[3]),max(bounds[4],extrema_now[4])]
            all(isfinite,extrema_now) && extrema_now[1] > 0 && 0 <= extrema_now[3] <= extrema_now[4] < 1 || error("Invalid state bounds")
            if f.integrator.time % 3600 == 0 || step == steps
                println(hourly,join(vcat(f.integrator.time,extrema_now,final[1:2],rain_current,
                    rates["bm_rain"],rates["lscale_rain"]),'\t'))
                flush(hourly)
            end
            if f.integrator.time % 21600 == 0 || step == steps
                @printf("dt=%d day=%.2f T=[%.2f,%.2f] qmax=%.6f BM/LC rain=%.4f/%.4f mm elapsed=%.1fs\n",
                    dt,f.integrator.time/86400,extrema_now[1],extrema_now[2],extrema_now[4],rain_current[1],rain_current[2],(time_ns()-start)/1.0e9)
                flush(stdout)
            end
        end
    finally
        close(hourly)
    end
    serialize(joinpath(output,"final.bin"),(u=copy(f.dyn.grid_u_c),v=copy(f.dyn.grid_v_c),t=copy(f.dyn.grid_t_c),
        q=copy(f.dyn.grid_q_c),ps=copy(f.dyn.grid_ps_c),weights=copy(f.mesh.wts),ph=copy(f.dyn.grid_p_half)))
    source_root = normpath(joinpath(@__DIR__,"..",".."))
    hashes = Dict(relpath(joinpath(folder,name),source_root) => bytes2hex(sha256(read(joinpath(folder,name))))
        for (folder,_,names) in walkdir(joinpath(source_root,"src")) for name in names if endswith(name,".jl"))
    hashes["exp/BettsMiller_Validation/validate.jl"] = bytes2hex(sha256(read(@__FILE__)))
    report = Dict("julia_version"=>string(VERSION),"threads"=>Threads.nthreads(),"days"=>days,"dt"=>dt,
        "steps"=>steps,"initial_state_sha256"=>initial_hash,"initial_ledger"=>initial,"final_ledger"=>dyn_ledger(f,:c),
        "bounds"=>bounds,"statistics"=>stats,"cumulative_process_changes"=>current,
        "accumulated_bm_rain_mm"=>rain_current[1],"accumulated_lscale_rain_mm"=>rain_current[2],
        "accumulated_evaporation_mm"=>rain_current[3],"source_sha256"=>hashes,
        "elapsed_seconds"=>(time_ns()-start)/1.0e9,"passed"=>true)
    open(joinpath(output,"result.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    return report
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 1 || error("Usage: validate.jl OUTPUT [DAYS=5] [DT=600]")
    run_validation(abspath(ARGS[1]),length(ARGS)>1 ? parse(Float64,ARGS[2]) : 5.0,
        length(ARGS)>2 ? parse(Int,ARGS[3]) : 600)
end
