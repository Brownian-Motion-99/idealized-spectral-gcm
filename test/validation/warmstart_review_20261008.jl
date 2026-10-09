# julia --project=. --startup-file=no --threads=1 \
#   test/validation/warmstart_review_20261008.jl NEW_OUTPUT [CHECKPOINT|none]
# Uses a new output directory; the supplied checkpoint is read only.
using Test, JGCM, JLD2, NCDatasets, TOML, SHA, Dates
include("../support/post_physics_state.jl")
const W = PostPhysicsStateSupport
const REVIEW_OUTPUT = abspath(get(ARGS, 1,
    "/tmp/jgcm_warmstart_" * Dates.format(now(), "yyyymmdd_HHMMSS")))
const REVIEW_CHECKPOINT = get(ARGS, 2,
    "/data92/garywu/undergrad_proposal/ctrl_BM/restart/restart_t315360000.jld2")
ispath(REVIEW_OUTPUT) && error("Use a new output directory: $REVIEW_OUTPUT")
mkpath(REVIEW_OUTPUT)
const FIELD_REPORT = open(joinpath(REVIEW_OUTPUT, "fields.tsv"), "w")
println(FIELD_REPORT, "phase\tfield\texpected_shape\tactual_shape\ttype_match\tbitwise_equal\tmax_abs_error\ttolerance\tfinite\tpassed")
const RESULTS = Dict{String,Any}("julia_version" => string(VERSION),
    "threads" => Threads.nthreads(), "output" => REVIEW_OUTPUT)
const STEP_TRACE = Dict{String,Vector{NamedTuple}}()
const FIRST_STEP_CHECK = Dict{String,Function}()

# Observe the actual driver before each step in this test process. Delegate to
# its original method, preserving the production time loop and dynamics.
function JGCM.Driver.Step_Dynamics!(config::Model_Config, mesh, atmo, dyn,
                                   integrator, semi, vert, physics_params)
    if haskey(STEP_TRACE, config.name)
        trace = STEP_TRACE[config.name]
        if isempty(trace) && haskey(FIRST_STEP_CHECK, config.name)
            FIRST_STEP_CHECK[config.name](dyn, integrator, semi)
        end
        push!(trace, (;time=integrator.time, start_time=integrator.start_time,
            end_time=integrator.end_time, dt=integrator.Δt,
            effective_dt=W.TI.Get_Δt(integrator), init_step=integrator.init_step))
    end
    invoke(JGCM.Driver.Step_Dynamics!, Tuple{Any,Any,Any,Any,Any,Any,Any,Any},
        config, mesh, atmo, dyn, integrator, semi, vert, physics_params)
end

array_fields(dyn) = filter(n -> getfield(dyn, n) isa Array, fieldnames(Dyn_Data))
bitwise_equal(a, b) = typeof(a) == typeof(b) && size(a) == size(b) &&
    reinterpret(UInt8, vec(a)) == reinterpret(UInt8, vec(b))

function compare_fields(phase, expected, actual; exact=false, reconstructed=false)
    same_bits = 0
    different_bits = String[]
    for name in array_fields(expected)
        a, b = getfield(expected, name), getfield(actual, name)
        shape_ok, type_ok = size(a) == size(b), typeof(a) == typeof(b)
        finite = all(isfinite, a) && all(isfinite, b)
        identical = bitwise_equal(a, b)
        err = shape_ok ? maximum(abs, a .- b) : Inf
        field_exact = exact || (reconstructed && name ∉ (:grid_vor, :grid_div))
        limit = field_exact ? 0.0 : 1e-18 + 5e-12 * maximum(abs, a)
        passed = shape_ok && type_ok && finite && (field_exact ? identical : err <= limit)
        @testset "$phase/$name" begin
            @test shape_ok
            @test type_ok
            @test finite
            @test passed
        end
        println(FIELD_REPORT, join((phase, name, size(a), size(b), type_ok,
            identical, err, limit, finite, passed), '\t'))
        identical ? (same_bits += 1) : push!(different_bits, string(name))
    end
    flush(FIELD_REPORT)
    Dict("fields" => length(array_fields(expected)), "bitwise_equal" => same_bits,
         "different_bits" => different_bits)
end

function read_state(path)
    dims, saved_time = jldopen(path, "r") do file
        @test file["restart_format_version"] == 3
        file["dimensions"], file["saved_time"]
    end
    dyn = Dyn_Data("review", dims...)
    @test Load_Restart_File!(dyn, path) === saved_time
    dyn, saved_time
end

function run_config(case, dir; nf=5, nlat=12, nd=8, steps=6, kwargs...)
    base = W.configuration(case; output=dir, nf, nlat, nd, steps)
    W.with_config(base; saving_frequency=600, output_interval=600,
        vars_to_output=[:u, :v, :t, :q, :ps, :bm_dt, :bm_dq, :bm_precip], kwargs...)
end

function run_driver(config)
    STEP_TRACE[config.name] = NamedTuple[]
    JGCM_Simulate(config)
    STEP_TRACE[config.name]
end

function check_clock(config, trace, start_time; warm)
    @test length(trace) == div(config.end_time, config.Δt)
    for (i, state) in enumerate(trace)
        @test state.time == start_time + (i - 1) * config.Δt
        @test state.start_time == start_time
        @test state.end_time == start_time + config.end_time
        @test state.dt == config.Δt
        @test state.init_step == (!warm && i == 1)
        @test state.effective_dt == ((!warm && i == 1) ? config.Δt : 2config.Δt)
    end
    # NetCDF records are labeled with the absolute time at interval end.
    all_times = Float64[]
    for i in eachindex(trace)
        chunk_start = start_time + (i - 1) * config.Δt
        chunk_path = splitext(config.output_filename)[1] * "_t$chunk_start.nc"
        NCDataset(chunk_path) do ds
            times = Float64.(ds["time"].var[:]) * config.day_to_sec
            @test length(times) == 1
            @test times[1] ≈ chunk_start + config.Δt atol=1e-6 rtol=0
            append!(all_times, times)
        end
    end
    @test length(all_times) == length(trace)
    [Dict(string(k) => v for (k, v) in pairs(state)) for state in trace]
end

function split_run(case; nf=5, nlat=12, nd=8, label=string(case))
    @testset "$label uninterrupted versus warm start" begin
        full_dir = joinpath(REVIEW_OUTPUT, label, "uninterrupted")
        cold_dir = joinpath(REVIEW_OUTPUT, label, "cold")
        warm_dir = joinpath(REVIEW_OUTPUT, label, "warm")
        full = run_config(case, full_dir; nf, nlat, nd, steps=6, name=label * "_full")
        cold = run_config(case, cold_dir; nf, nlat, nd, steps=3, name=label * "_cold")
        full_trace, cold_trace = run_driver(full), run_driver(cold)
        checkpoint = joinpath(cold_dir, "restart", "restart_t1800.jld2")
        source, source_time = read_state(checkpoint)
        @test source_time === Int64(1800)
        reference, _ = read_state(joinpath(full_dir, "restart", "restart_t1800.jld2"))
        info = Dict{String,Any}("dimensions" => [nf, nf+1, 2nlat, nlat, nd])
        info["cold_checkpoint"] = compare_fields(label * "_cold_checkpoint", reference, source)
        warm = run_config(case, warm_dir; nf, nlat, nd, steps=3,
            name=label * "_warm", is_restart=true, restart_file=checkpoint,
            initial_condition=(args...) -> error("Warm start called cold initialization"))
        FIRST_STEP_CHECK[warm.name] = (dyn, integrator, semi) -> begin
            # The only reconstructed fields are grid vorticity and divergence.
            expected = deepcopy(source)
            mesh = Spectral_Spherical_Mesh(nf, nf+1, 2nlat, nlat, nd, warm.radius)
            Trans_Spherical_To_Grid!(mesh, expected.spe_vor_c, expected.grid_vor)
            Trans_Spherical_To_Grid!(mesh, expected.spe_div_c, expected.grid_div)
            info["first_loaded_state"] = compare_fields(label * "_loaded", expected, dyn;
                reconstructed=true)
            @test integrator.time === Int64(1800)
            @test !integrator.init_step
        end
        warm_trace = run_driver(warm)
        info["full_clock"] = check_clock(full, full_trace, 0; warm=false)
        info["cold_clock"] = check_clock(cold, cold_trace, 0; warm=false)
        info["warm_clock"] = check_clock(warm, warm_trace, 1800; warm=true)
        for time in (2400, 3000, 3600)
            expected, expected_time = read_state(joinpath(full_dir, "restart", "restart_t$time.jld2"))
            actual, actual_time = read_state(joinpath(warm_dir, "restart", "restart_t$time.jld2"))
            @test expected_time === actual_time === Int64(time)
            info["t$time"] = compare_fields(label * "_t$time", expected, actual)
        end
        RESULTS[label] = info
    end
end

function actual_checkpoint_run(path)
    @testset "Configured T42 checkpoint" begin
        initial_hash = open(sha256, path)
        dims, source_time = jldopen(path, "r") do file
            file["dimensions"], file["saved_time"]
        end
        nf, ns, nlon, nlat, nd = dims
        source, loaded = Dyn_Data("source", dims...), Dyn_Data("loaded", dims...)
        buffers = Dict(n => getfield(loaded, n) for n in array_fields(loaded))
        for a in values(buffers); fill!(a, NaN); end
        jldopen(path, "r") do file
            @test Set(keys(file["state"])) == Set(string.(array_fields(source)))
            for n in array_fields(source)
                a = file["state/$n"]
                @test size(a) == size(getfield(source, n))
                @test typeof(a) == typeof(getfield(source, n))
                copyto!(getfield(source, n), a)
            end
        end
        @test Load_Restart_File!(loaded, path) === source_time
        @test source_time % 600 == 0
        @test all(getfield(loaded, n) === a for (n, a) in buffers)
        info = Dict{String,Any}("path" => path, "sha256" => bytes2hex(initial_hash),
            "dimensions" => collect(dims), "saved_time_seconds" => source_time,
            "dt_seconds" => 600, "warm_steps" => 4)
        info["raw_load"] = compare_fields("actual_raw_load", source, loaded; exact=true)

        # Independent in-memory continuation supplies a comparison after every
        # resumed step; the production driver performs the second continuation.
        reference_config = run_config(:bm_on, joinpath(REVIEW_OUTPUT, "actual_reference");
            nf, nlat, nd, steps=4, name="actual_reference", is_restart=true,
            restart_file=path, physics_params=W.parameters(:bm_on))
        reference = W.frame(reference_config; manufactured=false)
        reference.integrator.time = source_time
        reference.integrator.start_time = source_time
        reference.integrator.end_time = source_time + reference_config.end_time
        W.SI.Update_Init_Step!(reference.semi)
        info["initialized_load"] = compare_fields("actual_initialized_load", source,
            reference.dyn; reconstructed=true)
        warm_config = W.with_config(reference_config; name="actual_driver",
            physics_params=W.parameters(:bm_on),
            output_path=joinpath(REVIEW_OUTPUT, "actual_driver"),
            output_filename=joinpath(REVIEW_OUTPUT, "actual_driver", "output.nc"),
            logger=joinpath(REVIEW_OUTPUT, "actual_driver", "progress.log"),
            initial_condition=(args...) -> error("Warm start called cold initialization"))
        FIRST_STEP_CHECK[warm_config.name] = (dyn, integrator, semi) -> begin
            info["driver_loaded_state"] = compare_fields("actual_driver_loaded",
                reference.dyn, dyn; reconstructed=true)
            @test semi.wave_matrix ≈ reference.semi.wave_matrix rtol=5e-13
            @test W.TI.Get_Δt(integrator) == 1200
        end
        trace = run_driver(warm_config)
        info["clock"] = check_clock(warm_config, trace, source_time; warm=true)
        for step in 1:4
            W.step!(reference)
            time = source_time + 600step
            actual, saved_time = read_state(joinpath(warm_config.output_path,
                "restart", "restart_t$time.jld2"))
            @test saved_time === Int64(time)
            info["step$step"] = compare_fields("actual_step$step", reference.dyn, actual)
        end
        @test open(sha256, path) == initial_hash
        info["source_checkpoint_unchanged"] = true
        RESULTS["actual_checkpoint"] = info
    end
end

try
    @testset "Warm start review" begin
        include("../test_Memory_And_Restart.jl")
        split_run(:bm_off)
        split_run(:bm_on)
        split_run(:bm_on; nf=42, nlat=64, nd=20, label="t42_bm_on")
        REVIEW_CHECKPOINT == "none" || actual_checkpoint_run(REVIEW_CHECKPOINT)
    end
    RESULTS["passed"] = true
finally
    close(FIELD_REPORT)
    root = normpath(joinpath(@__DIR__, "../.."))
    source_files = sort([relpath(joinpath(dir, file), root)
        for (dir, _, files) in walkdir(joinpath(root, "src"))
        for file in files if endswith(file, ".jl")])
    append!(source_files, ["test/support/post_physics_state.jl",
        "test/test_Memory_And_Restart.jl", "test/validation/warmstart_review_20261008.jl"])
    RESULTS["source_sha256"] = Dict(file => bytes2hex(open(sha256, joinpath(root, file)))
        for file in source_files)
    open(joinpath(REVIEW_OUTPUT, "results.toml"), "w") do io
        TOML.print(io, RESULTS; sorted=true)
    end
end
println("Review artifacts: ", REVIEW_OUTPUT)
