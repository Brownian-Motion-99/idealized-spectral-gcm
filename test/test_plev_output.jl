"""
test_plev_output.jl

Quick test: PrimitiveEquation with do_plev_output=true.
Uses Held-Suarez forcing (always numerically stable) at low resolution for speed.
Verifies that both ds_raw (output_t0.nc) and ds_plev (output_t0_plev.nc) are created.
"""

using JGCM
using Test

experiment_name  = "PE_plev_test"
output_path_base = joinpath("exp", "test_tmp", experiment_name)
isdir(output_path_base) && rm(output_path_base; recursive=true)
mkpath(output_path_base)

config = Model_Config(
    name        = experiment_name,
    model_type  = :PrimitiveEquation,
    institution = "Test Suite",

    # T21 with 10 levels — small enough to run fast, large enough to be stable
    num_fourier = 21,
    nθ          = 32,
    nd          = 10,

    radius     = 6371.0e3,
    omega      = 7.292e-5,
    grav       = 9.80,
    day_to_sec = 86400,

    vert_coord_option      = "even_sigma",
    vert_difference_option = "simmons_and_burridge",
    vert_ref_level_option  = "second_centered_wts",

    Δt       = 600,
    end_time = 86400,       # 1-day run

    damping_order = 4,
    damping_coef  = 1.15741e-4,
    robert_coef   = 0.04,
    implicit_coef = 0.5,

    is_restart        = false,
    restart_file      = "",
    saving_frequency  = 0,   # no rotation for this test

    moisture_processes = false,

    output_path     = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger          = joinpath(output_path_base, "logger.log"),

    vars_to_output  = [:u, :v, :t, :ps],
    output_interval = 3600,

    do_plev_output  = true,
    pressure_levels = [85000.0, 50000.0, 20000.0],

    initial_condition = :Moist_Spinup,
    physics_params    = Dict{String, Any}(
        "do_mass_correction"      => true,
        "do_energy_correction"    => true,
        "do_water_correction"     => false,
        "use_virtual_temperature" => false,
        # Held-Suarez forcing keeps the model stable
        "do_HS_Forcing" => true,
        "σ_b"           => 0.7,
        "k_a"           => 1.0 / 40.0,
        "k_s"           => 1.0 / 4.0,
        "k_f"           => 1.0 / 1.0,
        "ΔT_y"          => 60.0,
        "Δθ_z"          => 10.0,
        # Disable moist physics
        "do_Lscale_Cond"          => false,
        "do_Sensible_Heating"     => false,
        "do_Surface_Evaporation"  => false,
        "do_Implicit_PBL_Scheme"  => false,
    )
)

@testset "PrimitiveEquation with do_plev_output=true" begin
    JGCM_Simulate(config)

    base = splitext(config.output_filename)[1]

    @info "Files present: $(readdir(output_path_base))"

    # Native grid output (always)
    @test isfile("$(base)_t0.nc")

    # Pressure-level output (optional, requested here)
    @test isfile("$(base)_t0_plev.nc")

    # No plain output.nc (new naming convention)
    @test !isfile("$(base).nc")
end

@info "PrimitiveEquation plev output test complete."
