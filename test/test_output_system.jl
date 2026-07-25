"""
test_output_system.jl

End-to-end tests for the chunked NC output + saving_frequency system.

Test 1  — Barotropic, 2-day run with saving_frequency=86400
          Expects: output_t0.nc  (chunk 0)
                   output_t86400.nc  (chunk 1, created at rotation)
                   restart/restart_t86400.jld2  (JLD2 checkpoint)
                   do_plev_output is ignored for 2D → no _plev files

Test 2  — Sanity-check: do_plev_output=true with pressure_levels=[] must throw an error.

Test 3  — saving_frequency=0 (disabled): single output_t0.nc produced, no rotation.
"""

using JGCM
using Test
using NCDatasets

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function make_barotropic_config(;
    experiment_name,
    end_days,
    saving_frequency,
    do_plev_output = false,
    pressure_levels = Float64[],
    output_interval = 3600
)
    output_path_base = joinpath("exp", "test_tmp", experiment_name)
    mkpath(output_path_base)

    Model_Config(
        name        = experiment_name,
        model_type  = :Barotropic,
        institution = "Test Suite",

        # Low-resolution T21 for speed
        num_fourier = 21,
        nθ          = 32,
        nd          = 1,

        radius     = 6371.0e3,
        omega      = 7.292e-5,
        grav       = 9.80,
        day_to_sec = 86400,

        vert_coord_option      = "even_sigma",
        vert_difference_option = "none",
        vert_ref_level_option  = "none",

        Δt       = 1200,
        end_time = 86400 * end_days,

        damping_order = 4,
        damping_coef  = 1.0e-4,
        robert_coef   = 0.04,
        implicit_coef = 0.0,

        is_restart        = false,
        restart_file      = "",
        saving_frequency  = saving_frequency,

        moisture_processes = false,
        num_tracers        = 0,

        output_path     = output_path_base,
        output_filename = joinpath(output_path_base, "output.nc"),
        logger          = joinpath(output_path_base, "logger.log"),

        vars_to_output  = [:u, :v, :vor],
        output_interval = output_interval,

        do_plev_output  = do_plev_output,
        pressure_levels = pressure_levels,

        initial_condition = :Barotropic_Jet,
        physics_params    = Dict{String, Any}()
    )
end

# ---------------------------------------------------------------------------
# Test 1: Chunked output — 2-day Barotropic run, rotate at day 1
# ---------------------------------------------------------------------------
@testset "Test 1: Barotropic chunked NC output (saving_frequency=86400)" begin
    config = make_barotropic_config(
        experiment_name  = "barotropic_chunk_test",
        end_days         = 2,
        saving_frequency = 86400,
        output_interval  = 3600
    )

    base = splitext(config.output_filename)[1]    # strip ".nc"
    restart_dir = joinpath(config.output_path, "restart")

    # Clean up any previous test artifacts
    isdir(config.output_path) && rm(config.output_path; recursive=true)
    mkpath(config.output_path)

    JGCM_Simulate(config)

    @info "Checking output files in: $(config.output_path)"
    @info "Files present: $(readdir(config.output_path))"

    @test isfile("$(base)_t0.nc")
    @test isfile("$(base)_t86400.nc")
    @test !isfile("$(base).nc")

    # JLD2 checkpoint should also exist
    @test isfile(joinpath(restart_dir, "restart_t86400.jld2"))

    # No plev files expected for Barotropic
    @test !isfile("$(base)_t0_plev.nc")
end

# ---------------------------------------------------------------------------
# Test 2: Sanity check — do_plev_output=true with empty pressure_levels
# ---------------------------------------------------------------------------
@testset "Test 2: Sanity check — do_plev_output without pressure_levels throws" begin
    experiment_name  = "sanity_check_plev"
    output_path_base = joinpath("exp", "test_tmp", experiment_name)
    isdir(output_path_base) && rm(output_path_base; recursive=true)
    mkpath(output_path_base)

    config_bad = Model_Config(
        name        = experiment_name,
        model_type  = :PrimitiveEquation,
        institution = "Test Suite",

        num_fourier = 10,
        nθ          = 16,
        nd          = 5,

        radius     = 6371.0e3,
        omega      = 7.292e-5,
        grav       = 9.80,
        day_to_sec = 86400,

        vert_coord_option      = "even_sigma",
        vert_difference_option = "simmons_and_burridge",
        vert_ref_level_option  = "second_centered_wts",

        Δt       = 1200,
        end_time = 86400,

        damping_order = 4,
        damping_coef  = 1.15741e-4,
        robert_coef   = 0.04,
        implicit_coef = 0.5,

        is_restart        = false,
        restart_file      = "",
        saving_frequency  = 0,

        moisture_processes = true,
        num_tracers        = 1,

        output_path     = output_path_base,
        output_filename = joinpath(output_path_base, "output.nc"),
        logger          = joinpath(output_path_base, "logger.log"),

        vars_to_output  = [:u, :v, :t, :ps],
        output_interval = 3600,

        # The problematic combination: request plev output but forget to set levels
        do_plev_output  = true,
        pressure_levels = Float64[],   # empty → should throw

        initial_condition = :Moist_Spinup,
        physics_params    = Dict{String, Any}(
            "do_mass_correction"      => false,
            "do_energy_correction"    => false,
            "do_water_correction"     => false,
            "use_virtual_temperature" => false,
            "do_Lscale_Cond"          => false,
            "do_Sensible_Heating"     => false,
            "do_Surface_Evaporation"  => false,
            "do_Implicit_PBL_Scheme"  => false,
            "do_HS_Forcing"           => false,
        )
    )

    @test_throws ErrorException JGCM_Simulate(config_bad)
end

# ---------------------------------------------------------------------------
# Test 3: saving_frequency=0 — no rotation, single chunk per run
# ---------------------------------------------------------------------------
@testset "Test 3: Barotropic single-chunk run (saving_frequency=0)" begin
    config = make_barotropic_config(
        experiment_name  = "barotropic_nochunk_test",
        end_days         = 1,
        saving_frequency = 0,
        output_interval  = 3600
    )

    base = splitext(config.output_filename)[1]
    restart_dir = joinpath(config.output_path, "restart")

    isdir(config.output_path) && rm(config.output_path; recursive=true)
    mkpath(config.output_path)

    JGCM_Simulate(config)

    @info "Checking output files in: $(config.output_path)"
    @info "Files present: $(readdir(config.output_path))"

    @test isfile("$(base)_t0.nc")
    @test !isfile("$(base)_t86400.nc")

    # No JLD2 checkpoint when saving_frequency=0
    no_jld2 = !isdir(restart_dir) || isempty(filter(f -> endswith(f, ".jld2"), readdir(restart_dir)))
    @test no_jld2
end

# ---------------------------------------------------------------------------
# Test 4: Warm-start with run duration shorter than saving_frequency
#
# Phase 1 (cold start): 2-day run, saving_frequency=86400 (1 day).
#   → Produces restart_t86400.jld2 and restart_t172800.jld2.
#   → NC chunks: output_t0.nc (day 0–1) and output_t86400.nc (day 1–2).
#
# Phase 2 (warm start): restart from t=172800 (day 2), end_time=86400 (1 day),
#   saving_frequency=172800 (2 days) — next checkpoint would be at t=345600,
#   which is beyond the segment end (t=259200). No rotation fires.
#   → Single NC chunk: output_t172800.nc covering day 2–3 (24 time points).
#   → No new JLD2 checkpoint written.
# ---------------------------------------------------------------------------
@testset "Test 4: Warm-start output — duration shorter than saving_frequency" begin
    cold_path = joinpath("exp", "test_tmp", "warmstart_cold_phase")
    warm_path = joinpath("exp", "test_tmp", "warmstart_warm_phase")

    isdir(cold_path) && rm(cold_path; recursive=true)
    isdir(warm_path) && rm(warm_path; recursive=true)
    mkpath(cold_path)
    mkpath(warm_path)

    # ------------------------------------------------------------------
    # Phase 1: cold start, 2 days, checkpoint every 1 day
    # ------------------------------------------------------------------
    config_cold = make_barotropic_config(
        experiment_name  = "warmstart_cold_phase",
        end_days         = 2,
        saving_frequency = 86400,   # checkpoint + rotation every 1 day
        output_interval  = 3600     # 24 snapshots per day
    )

    JGCM_Simulate(config_cold)

    cold_base    = splitext(config_cold.output_filename)[1]
    restart_dir  = joinpath(cold_path, "restart")
    restart_file = joinpath(restart_dir, "restart_t172800.jld2")

    # Phase 1 sanity checks
    @test isfile("$(cold_base)_t0.nc")
    @test isfile("$(cold_base)_t86400.nc")
    @test isfile(joinpath(restart_dir, "restart_t86400.jld2"))
    @test isfile(restart_file)

    # ------------------------------------------------------------------
    # Phase 2: warm start from day 2, run for 1 day
    #   saving_frequency=172800 (2 days) → next hit at t=345600 > t=259200
    #   so no chunk rotation and no checkpoint should be produced
    # ------------------------------------------------------------------
    config_warm = Model_Config(
        name        = "warmstart_warm_phase",
        model_type  = :Barotropic,
        institution = "Test Suite",

        num_fourier = 21,
        nθ          = 32,
        nd          = 1,

        radius     = 6371.0e3,
        omega      = 7.292e-5,
        grav       = 9.80,
        day_to_sec = 86400,

        vert_coord_option      = "even_sigma",
        vert_difference_option = "none",
        vert_ref_level_option  = "none",

        Δt       = 1200,
        end_time = 86400,           # 1-day segment

        damping_order = 4,
        damping_coef  = 1.0e-4,
        robert_coef   = 0.04,
        implicit_coef = 0.0,

        is_restart       = true,
        restart_file     = restart_file,
        saving_frequency = 172800,  # 2 days — will not trigger during 1-day run

        moisture_processes = false,
        num_tracers        = 0,

        output_path     = warm_path,
        output_filename = joinpath(warm_path, "output.nc"),
        logger          = joinpath(warm_path, "logger.log"),

        vars_to_output  = [:u, :v, :vor],
        output_interval = 3600,

        do_plev_output  = false,
        pressure_levels = Float64[],

        initial_condition = :Barotropic_Jet,
        physics_params    = Dict{String, Any}()
    )

    JGCM_Simulate(config_warm)

    warm_base        = splitext(config_warm.output_filename)[1]
    warm_restart_dir = joinpath(warm_path, "restart")

    # Single NC chunk starting at restart time (t=172800 s = day 2)
    @test isfile("$(warm_base)_t172800.nc")

    # No rotation fired — no chunk at the segment end
    @test !isfile("$(warm_base)_t259200.nc")

    # No checkpoint written during the 1-day warm-start run
    no_new_jld2 = !isdir(warm_restart_dir) ||
                  isempty(filter(f -> endswith(f, ".jld2"), readdir(warm_restart_dir)))
    @test no_new_jld2

    # NC content: 24 time points covering day 2 → day 3
    # Use .var[:] to get raw Float64 values (bypasses DateTime calendar decoding)
    NCDatasets.Dataset("$(warm_base)_t172800.nc", "r") do ds
        t = ds["time"].var[:]
        @test length(t) == 24                          # 86400 s / 3600 s per snapshot
        @test t[1]   ≈ (172800 + 3600) / 86400        # first flush at day 2 + 1 h
        @test t[end] ≈ 259200 / 86400                 # last flush at day 3
    end
end

@info "All output system tests complete."
