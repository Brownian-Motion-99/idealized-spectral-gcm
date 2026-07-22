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

@info "All output system tests complete."
