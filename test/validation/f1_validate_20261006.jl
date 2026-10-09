# Run in either the preserved baseline or repaired source tree:
# julia --project=. --compiled-modules=no --threads=1 \
#   test/validation/f1_validate_20261006.jl OUTPUT small|t42
# OUTPUT must be a new directory. Generated binary states stay outside the repo.
using JGCM, Serialization, TOML, SHA, Statistics, Printf
include("../support/post_physics_state.jl")
const S = PostPhysicsStateSupport
const ROOT = normpath(joinpath(@__DIR__,"../.."))

function ledger(f)
    d, m, a, v = f.dyn, f.mesh, f.atmo, f.vert
    water, energy, dry = 0.0, 0.0, 0.0
    for k in 1:m.nd, j in 1:m.nθ, i in 1:m.nλ
        mass = (v.Δak[k]+v.Δbk[k]*d.grid_ps_c[i,j,1])*m.wts[j]/(2m.nλ*a.grav)
        q = d.grid_q_c[i,j,k]
        water += mass*q
        dry += mass*(1-q)
        energy += mass*(a.cp_air*d.grid_t_c[i,j,k]+a.Lv*q+
            0.5*(d.grid_u_c[i,j,k]^2+d.grid_v_c[i,j,k]^2))
    end
    (;water,energy,dry)
end

function run_case(output,case; nf=5,nlat=12,nd=8,steps=2,
                  manufactured=true,manual_refresh=false)
    config = S.configuration(case;output,nf,nlat,nd,steps)
    effective = Dict(string(key)=>getfield(config,key) for key in fieldnames(Model_Config))
    effective["model_type"] = string(config.model_type)
    effective["initial_condition"] = string(config.initial_condition)
    effective["vars_to_output"] = string.(config.vars_to_output)
    effective["physics_params"] = copy(config.physics_params)
    effective["physics_params"]["PBL_Top_Mode"] = string(config.physics_params["PBL_Top_Mode"])
    f = S.frame(config;manufactured)
    vor,div = similar(f.dyn.grid_vor),similar(f.dyn.grid_div)
    initial_ledger = ledger(f)
    stats = Dict{String,Any}(
        "effective_configuration"=>effective,"manufactured_winds"=>manufactured,
        "manual_refresh"=>manual_refresh,"steps"=>steps,
        "max_vor_error"=>0.0,"max_div_error"=>0.0,
        "max_vor_relative_error"=>0.0,"max_div_relative_error"=>0.0,
        "invariant_failures"=>0,"bm_active_steps"=>0,
        "max_wind"=>0.0,"min_temperature"=>Inf,"max_temperature"=>-Inf,
        "min_pressure"=>Inf,"min_humidity"=>Inf,"max_humidity"=>-Inf,
        "bm_rain_rate_integral"=>0.0,"total_rain_rate_integral"=>0.0,
        "initial_water"=>initial_ledger.water,"initial_energy"=>initial_ledger.energy,
        "initial_dry_mass"=>initial_ledger.dry,
    )
    times = Float64[]
    tag = string(case)*(manual_refresh ? "_reference" : "")
    observation = open(joinpath(output,tag*"_steps.tsv"),"w")
    println(observation,"step\ttime_s\tvor_error_s-1\tdiv_error_s-1\twater_kg_m-2\tenergy_J_m-2\tmean_ps_Pa\tbm_rain_kg_m-2_s-1\tstep_seconds")
    try
        for step in 1:steps
            elapsed = @elapsed S.step!(f)
            if manual_refresh
                Trans_Spherical_To_Grid!(f.mesh,f.dyn.spe_vor_c,f.dyn.grid_vor)
                Trans_Spherical_To_Grid!(f.mesh,f.dyn.spe_div_c,f.dyn.grid_div)
            end
            if step <= 2
                serialize(joinpath(output,"$(tag)_step$(step).bin"),S.snapshot(f.dyn))
            end
            e = S.errors!(f,vor,div)
            stats["max_vor_error"] = max(stats["max_vor_error"],e.evor)
            stats["max_div_error"] = max(stats["max_div_error"],e.ediv)
            stats["max_vor_relative_error"] = max(stats["max_vor_relative_error"],e.evor/max(e.svor,1e-30))
            stats["max_div_relative_error"] = max(stats["max_div_relative_error"],e.ediv/max(e.sdiv,1e-30))
            if e.evor > 1e-18+5e-13*e.svor || e.ediv > 1e-18+5e-13*e.sdiv
                stats["invariant_failures"] += 1
            end
            d = f.dyn
            for key in (:grid_u_c,:grid_v_c,:grid_t_c,:grid_q_c,:grid_ps_c)
                all(isfinite,getfield(d,key)) || error("Nonfinite $key at step $step")
            end
            minimum(d.grid_t_c)>0 && minimum(d.grid_ps_c)>0 || error("Nonpositive T/ps")
            minimum(d.grid_q_c)>=-1e-14 && maximum(d.grid_q_c)<1 || error("Invalid q")
            stats["max_wind"] = max(stats["max_wind"],maximum(hypot.(d.grid_u_c,d.grid_v_c)))
            stats["min_temperature"] = min(stats["min_temperature"],minimum(d.grid_t_c))
            stats["max_temperature"] = max(stats["max_temperature"],maximum(d.grid_t_c))
            stats["min_pressure"] = min(stats["min_pressure"],minimum(d.grid_ps_c))
            stats["min_humidity"] = min(stats["min_humidity"],minimum(d.grid_q_c))
            stats["max_humidity"] = max(stats["max_humidity"],maximum(d.grid_q_c))
            any(!iszero,d.grid_bm_t_tendency) && (stats["bm_active_steps"] += 1)
            bm = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(f.mesh,d.grid_bm_precip)
            rain = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(f.mesh,d.grid_precip)
            stats["bm_rain_rate_integral"] += config.Δt*bm
            stats["total_rain_rate_integral"] += config.Δt*rain
            l = ledger(f)
            ps = JGCM.Spectral_Spherical_Mesh_Module.Area_Weighted_Global_Mean(f.mesh,d.grid_ps_c)
            println(observation,join((step,f.integrator.time,e.evor,e.ediv,l.water,l.energy,ps,bm,elapsed),'\t'))
            step > 4 && push!(times,elapsed)
        end
    finally
        close(observation)
    end
    stats["warmed_median_step_seconds"] = isempty(times) ? 0.0 : median(times)
    stats["warmed_mean_step_seconds"] = isempty(times) ? 0.0 : mean(times)
    final = ledger(f)
    stats["final_water"] = final.water; stats["final_energy"] = final.energy
    stats["final_dry_mass"] = final.dry
    # Warmed allocation/cost probe excludes the validation diagnostics.
    stats["warmed_step_allocations_bytes"] = @allocated S.step!(f)
    if isdefined(S.SD,:_refresh_grid_vor_div!)
        S.SD._refresh_grid_vor_div!(f.mesh,f.dyn.spe_vor_c,f.dyn.spe_div_c,
            f.dyn.grid_vor,f.dyn.grid_div)
        stats["refresh_allocations_bytes"] = @allocated S.SD._refresh_grid_vor_div!(
            f.mesh,f.dyn.spe_vor_c,f.dyn.spe_div_c,f.dyn.grid_vor,f.dyn.grid_div)
        samples = [@elapsed S.SD._refresh_grid_vor_div!(f.mesh,f.dyn.spe_vor_c,
            f.dyn.spe_div_c,f.dyn.grid_vor,f.dyn.grid_div) for _ in 1:20]
        stats["refresh_median_seconds"] = median(samples)
    end
    @printf("%s T%dL%d: %d steps, errors %.3e/%.3e, BM active %d, %.4f s/step\n",
        tag,nf,nd,steps,stats["max_vor_error"],stats["max_div_error"],
        stats["bm_active_steps"],stats["warmed_median_step_seconds"])
    stats
end

function main()
    length(ARGS)==2 || error("Usage: f1_validate_20261006.jl NEW_OUTPUT small|t42")
    output, mode = ARGS
    mode in ("small","t42") || error("Unknown mode $mode")
    ispath(output) && error("Output must be a new path: $output")
    mkpath(output)
    sources = Dict{String,String}()
    for (dir,_,files) in walkdir(joinpath(ROOT,"src")), name in files
        endswith(name,".jl") || continue
        path = joinpath(dir,name)
        sources[relpath(path,ROOT)] = bytes2hex(sha256(read(path)))
    end
    for name in ("exp/HSt42/HS.jl","test/support/post_physics_state.jl",
                 "test/validation/f1_validate_20261006.jl")
        sources[name] = bytes2hex(sha256(read(joinpath(ROOT,name))))
    end
    result = Dict{String,Any}("julia"=>string(VERSION),"threads"=>Threads.nthreads(),
        "source_root"=>ROOT,"sources"=>sources,"mode"=>mode)
    if mode=="small"
        result["hs"] = run_case(output,:hs)
        result["hs_reference"] = run_case(output,:hs;manual_refresh=true)
        result["none"] = run_case(output,:none)
    else
        for case in (:bm_off,:bm_on)
            result[string(case)] = run_case(output,case;
                nf=42,nlat=64,nd=20,steps=144,manufactured=false)
        end
    end
    open(joinpath(output,"results.toml"),"w") do io; TOML.print(io,result); end
end

main()
