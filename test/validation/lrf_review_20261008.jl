# Review evidence, not a replacement for the regression suite.
# julia --project=. --compiled-modules=no --threads=1 \
#   test/validation/lrf_review_20261008.jl contracts|online NEW_OUTPUT
# Artifacts and checkpoints are read only; generated files go to NEW_OUTPUT.
using Test, JGCM, JLD2, TOML, SHA, Statistics, Printf, NCDatasets
include("../support/post_physics_state.jl")
const S = PostPhysicsStateSupport
const ROOT = normpath(joinpath(@__DIR__, "../.."))
const ARTIFACT_ROOT = get(ENV, "JGCM_REVIEW_ARTIFACTS",
    "/home/garywu/undergrad_proposal/LRF/data")
const CONTROL_ROOT = get(ENV, "JGCM_REVIEW_CONTROLS",
    "/data92/garywu/undergrad_proposal")

function zero_alpha(state::Latitude_LRF_State)
    Latitude_LRF_State(state.B_by_lat, state.chi_reference, state.q0, 0.0, state.nλ)
end
function zero_alpha(state::Regularized_LRF_State)
    Regularized_LRF_State(state.B, state.chi_reference, state.taper,
        state.q0, 0.0, state.nλ)
end

function contracts(output)
    result = Dict{String,Any}()
    mesh = Spectral_Spherical_Mesh(42, 43, 128, 64, 20, 6.371e6)
    latitude = rad2deg.(mesh.θc)
    vert = Vert_Coordinate(128, 64, 20, "even_sigma",
        "simmons_and_burridge", "second_centered_wts")
    nominal_pfull = zeros(1, 1, 20)
    JGCM.Vertical_Interpolation_Module.Compute_Pressure_Grid!(
        nominal_pfull, vert.ak, vert.bk, fill(vert.p_ref, 1, 1))
    nominal_phalf = (vert.ak .+ vert.bk .* vert.p_ref) ./ 100.0
    for filename in ("latitude_lrf.jld2", "regularized_tapered_lrf.jld2")
        path = joinpath(ARTIFACT_ROOT, filename)
        data = JLD2.load(path)
        info = Dict{String,Any}("path" => path,
            "sha256" => bytes2hex(sha256(read(path))), "scheme" => data["scheme"],
            "alpha" => data["alpha"], "q0" => data["q0_kg_kg"],
            "source_sha256" => data["source_sha256"])
        @testset "$filename actual artifact" begin
            state = Load_LRF_State(path, 128, 64, 20; latitude)
            state1 = Load_LRF_State(path, 1, 64, 20; latitude)
            actual = similar(data["validation_q"])
            LRF!(state1, data["validation_q"], actual, 86400)
            err = maximum(abs, actual .- data["validation_tendency_K_s"])
            info["cross_language_max_error_K_s"] = err
            @test err <= 5e-14 * max(maximum(abs, data["validation_tendency_K_s"]), 1e-12)
            info["nominal_pfull_max_error_hPa"] =
                maximum(abs, data["pfull"] .- vec(nominal_pfull) ./ 100.0)
            info["nominal_phalf_max_error_hPa"] =
                maximum(abs, data["phalf"] .- nominal_phalf)
            @test info["nominal_pfull_max_error_hPa"] < 1e-9
            @test info["nominal_phalf_max_error_hPa"] < 1e-9
            reference = exp.(state.chi_reference) .- state.q0
            @test minimum(reference) >= 0.0
            q = repeat(permutedims(reshape(reference, 20, 64, 1), (3, 2, 1)), 128, 1, 1)
            tendency = similar(q)
            LRF!(state, q, tendency, 86400)
            info["reference_max_tendency_K_s"] = maximum(abs, tendency)
            @test maximum(abs, tendency) < 1e-17
            q .*= 1.05
            before = copy(q)
            LRF!(state, q, tendency, 86400)
            @test q == before
            LRF!(zero_alpha(state), q, tendency, 86400)
            @test all(iszero, tendency)
            fill!(q, 0.0)
            LRF!(state, q, tendency, 86400)
            @test all(isfinite, tendency)
            info["all_zero_humidity_tendency_range_K_day"] =
                [minimum(tendency)*86400, maximum(tendency)*86400]

            # Demonstrate an accepted mismatch, rather than certify this behavior.
            wrong = copy(data)
            wrong["pfull"] = reverse(data["pfull"])
            wrong["phalf"] = reverse(data["phalf"])
            wrong_path = joinpath(output, filename * ".wrong_pressure.jld2")
            jldopen(wrong_path, "w") do file
                for (key, value) in wrong; file[key] = value; end
            end
            accepted = Load_LRF_State(wrong_path, 128, 64, 20; latitude)
            @test typeof(accepted) == typeof(state)
            info["reversed_pressure_metadata_accepted"] = true

            # Isolate the actual artifact in the physics interface at both intervals.
            config = S.configuration(:none; nf=42, nlat=64, nd=20, steps=2)
            config.physics_params["do_LRF"] = true
            config.physics_params["LRF_state"] = state
            f = S.frame(config; manufactured=false)
            LRF!(state, f.dyn.grid_q_c, tendency, 86400)
            info["cold_start_tendency_range_K_day"] =
                [minimum(tendency)*86400, maximum(tendency)*86400]
            q .= repeat(permutedims(reshape(reference, 20, 64, 1), (3, 2, 1)), 128, 1, 1)
            q .*= 1.05
            LRF!(state, q, tendency, 86400)
            for initial in (true, false)
                f.integrator.init_step = initial
                fill!(f.dyn.grid_ps_n, 1e5)
                fill!(f.dyn.grid_u_n, 0.0); fill!(f.dyn.grid_v_n, 0.0)
                fill!(f.dyn.grid_t_n, 280.0); f.dyn.grid_q_n .= q
                S.PG.Pressure_Variables!(f.vert, f.dyn.grid_ps_n,
                    f.dyn.grid_p_half, f.dyn.grid_Δp, f.dyn.grid_lnp_half,
                    f.dyn.grid_p_full, f.dyn.grid_lnp_full)
                effective_dt = initial ? 600 : 1200
                JGCM.Atmos_Param_Module.Spectral_Physics!(config,
                    f.mesh, f.vert, f.atmo, f.dyn, f.semi, config.physics_params)
                @test maximum(abs, f.dyn.grid_t_n .-
                    (280.0 .+ effective_dt .* tendency)) <= 1e-12
                @test f.dyn.grid_lrf_tendency ≈ tendency rtol=1e-12
                @test f.dyn.grid_q_n ≈ q rtol=1e-12
                energy_before = S.SD.Compute_Corrections_Init(f.mesh, f.vert, f.atmo,
                    f.dyn.grid_ps_n, f.dyn.grid_energy_full, f.dyn.grid_u_n,
                    f.dyn.grid_v_n, f.dyn.grid_t_n, f.dyn.grid_q_n)[2]
                S.SD._synchronize_physics_next!(f.mesh, f.vert, f.atmo, f.dyn)
                energy_after = S.SD.Compute_Corrections_Init(f.mesh, f.vert, f.atmo,
                    f.dyn.grid_ps_n, f.dyn.grid_energy_full, f.dyn.grid_u_n,
                    f.dyn.grid_v_n, f.dyn.grid_t_n, f.dyn.grid_q_n)[2]
                @test energy_after ≈ energy_before rtol=1e-12
            end
            info["startup_and_leapfrog_temperature_increment_verified"] = true
            info["postphysics_energy_target_preserved"] = true
        end
        result[filename] = info
    end
    @testset "Legacy validation gaps demonstrated" begin
        state = LRF_State(fill(NaN, 2, 2, 2), zeros(2, 2, 2))
        tendency = zeros(2, 2, 2)
        LRF!(state, zeros(2, 2, 2), tendency, 86400)
        @test all(isnan, tendency)
        result["legacy_nonfinite_coefficients_accepted"] = true
        path = joinpath(output, "legacy_wrong_latitude.jld2")
        JLD2.jldsave(path; LRF_LW_q=ones(2, 2, 2), ref_q=zeros(2, 2, 2),
            latitude=[30.0, -30.0])
        @test Load_LRF_State(path, 2, 2, 2; latitude=[-30.0, 30.0]) isa LRF_State
        result["legacy_reversed_latitude_accepted"] = true
    end
    result
end

function output_stats(path)
    NCDataset(path) do ds
        info = Dict{String,Any}("path"=>path, "times_day"=>Float64.(ds["time"].var[:]),
            "global_attributes"=>Dict(string(k)=>string(v) for (k,v) in ds.attrib))
        for name in ("ua", "va", "ta", "hus", "ps", "lrf_dta_dt")
            variable = ds[name]
            values = Array(variable[ntuple(_ -> Colon(), ndims(variable))...])
            @test all(isfinite, values)
            info[name*"_range"] = [minimum(values), maximum(values)]
            if name == "lrf_dta_dt"
                # NCDatasets presents Julia axes (lon, lat, level, time).
                info["lrf_record_rms_K_day"] =
                    [sqrt(mean(abs2, values[:,:,:,i]))*86400 for i in axes(values,4)]
                info["lrf_record_max_abs_K_day"] =
                    [maximum(abs, values[:,:,:,i])*86400 for i in axes(values,4)]
            end
        end
        @test info["hus_range"][1] >= 0
        @test info["hus_range"][2] < 1
        @test info["ta_range"][1] > 100 && info["ta_range"][2] < 400
        @test info["ps_range"][1] > 0
        info
    end
end

function online(output)
    result = Dict{String,Any}()
    # Existing smoke configuration, full driver and NetCDF, with fresh output.
    ENV["JGCM_SMOKE_RESTART"] = joinpath(CONTROL_ROOT,"ctrl/restart/restart_t315360000.jld2")
    ENV["JGCM_SMOKE_OUTPUT"] = joinpath(output, "latitude_bm_off_5day")
    ENV["JGCM_SMOKE_LRF"] = "1"
    ENV["JGCM_LRF_FILE"] = joinpath(ARTIFACT_ROOT,"latitude_lrf.jld2")
    @testset "Five-day latitude-LRF driver smoke, BM off" begin
        path = get(ENV, "JGCM_REVIEW_BM_OFF_NC", "")
        if isempty(path)
            include(joinpath(ROOT, "exp/HSt42/LRF_latitude_smoke.jl"))
            path = joinpath(ENV["JGCM_SMOKE_OUTPUT"],"output_t315360000.nc")
        end
        result["latitude_bm_off_5day"] = output_stats(path)
    end
    # Current HSt42 physics parameters, from its own spun-up control checkpoint.
    for enabled in (false, true)
        tag = enabled ? "latitude_bm_on_1day" : "control_bm_on_1day"
        config = S.configuration(:bm_on; nf=42,nlat=64,nd=20,steps=144,
            output=joinpath(output,tag))
        config.physics_params["do_LRF"] = enabled
        config.physics_params["LRF_file"] = joinpath(ARTIFACT_ROOT,"latitude_lrf.jld2")
        config = S.with_config(config;
            is_restart=true,
            restart_file=joinpath(CONTROL_ROOT,"ctrl_BM/restart/restart_t315360000.jld2"),
            vars_to_output=[:u,:v,:t,:q,:ps,:lrf_dt,:bm_dt,:bm_dq,:bm_precip,:precip],
            output_interval=21600)
        @testset "$tag full driver" begin
            JGCM_Simulate(config)
            result[tag] = output_stats(joinpath(config.output_path,"output_t315360000.nc"))
        end
    end
    result
end

function main()
    length(ARGS)==2 || error("Usage: lrf_review_20261008.jl contracts|online NEW_OUTPUT")
    mode, output = ARGS
    mode in ("contracts","online") || error("Unknown mode $mode")
    ispath(output) && error("Output must be a new path")
    mkpath(output)
    sources = Dict{String,String}()
    for name in ("src/Physics/LRF.jl","src/Physics/Spectral_Physics_Interface.jl",
                 "src/Driver.jl","src/Dynamics/Spectral_Dynamics.jl",
                 "exp/HSt42/HS.jl","exp/HSt42/LRF_latitude_smoke.jl")
        sources[name] = bytes2hex(sha256(read(joinpath(ROOT,name))))
    end
    result = Dict{String,Any}("julia"=>string(VERSION),
        "threads"=>Threads.nthreads(),"mode"=>mode,"sources"=>sources)
    result["results"] = mode=="contracts" ? contracts(output) : online(output)
    open(joinpath(output,"results.toml"),"w") do io; TOML.print(io,result); end
end
main()
