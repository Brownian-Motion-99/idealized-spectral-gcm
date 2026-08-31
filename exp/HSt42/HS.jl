using JGCM

# Moist-physics switches
do_betts_miller = false
do_lscale_condensation = true

# LRF requires a resolution-matched JLD2 file. It is enabled when the file is
# supplied, for example:
#   JGCM_LRF_FILE=/path/to/LRF_state.jld2 julia --project exp/HSt42/HS.jl
lrf_file = get(ENV, "JGCM_LRF_FILE", "")
do_lrf = !isempty(lrf_file)

# 1. Physics Configuration
physics_params = Dict{String,Any}(

    # Corrections
    "do_mass_correction" => true,
    "do_energy_correction" => true,
    # Post-dynamics physics carries its own water budget; the grid tracer
    # correction preserves that post-physics integral.
    "do_water_correction" => true,
    "use_virtual_temperature" => true,

    # Betts-Miller convection
    "do_Betts_Miller" => do_betts_miller,
    "bm_tau" => 7200.0,
    "bm_relative_humidity" => 0.8,
    # Optional physical background humidity.
    "initial_humidity_floor" => 0.0,

    # Grid scale condensation
    "do_Lscale_Cond" => do_lscale_condensation,
    "condensation_heating_fraction" => 1.0,

    # Linear response function
    "do_LRF" => do_lrf,
    "LRF_file" => lrf_file,

    # PBL fluxes
    "do_Sensible_Heating" => true,
    "C_H" => 0.0044,
    "do_Surface_Evaporation" => true,
    "C_E" => 0.0044,
    "do_Implicit_PBL_Scheme" => true,
    "C_D" => 0.0044,

    # "PBL_Top_Mode"  => :ModelLevel,
    # "PBL_Top_Value" => 4,
    "PBL_Top_Mode" => :PressureLevel,
    "PBL_Top_Value" => 85000.0,

    # Held-Suarez
    "do_HS_Forcing" => true,
    "σ_b" => 0.7,
    "k_a" => 1.0 / (40.0),
    "k_s" => 1.0 / (4.0),
    "k_f" => 1.0 / (1.0),
    "T_equator" => 294.0,
    "T_stratosphere" => 200.0,
    "ΔT_y" => 65.0,
    "Δθ_z" => 10.0,

    # TODO: cumulus parameterization
    # TODO: radiation parameterization
)

# 2. Define Output Paths *Before* Configuration
experiment_name = "ctrl"
output_path_base = joinpath("/data92/garywu/undergrad_proposal", experiment_name)
mkpath(output_path_base)

# 3. Model Configuration
config = Model_Config(
    name = "HS_Moist_T42",
    institution = "Group of Chaos and Predictability, Department of Atmospheric Sciences, National Taiwan University",
    model_type = :PrimitiveEquation,

    # Resolution
    num_fourier = 42,
    nθ = 64,
    nd = 20,

    # Vertical Coordinate
    vert_coord_option = "even_sigma",
    vert_difference_option = "simmons_and_burridge",
    vert_ref_level_option = "second_centered_wts",

    # Planet Settings
    radius = 6371.0e3,
    omega = 7.292e-5,
    grav = 9.80,
    day_to_sec = 86400,

    # Time Integration
    Δt = 600,
    end_time = 86400 * 10000,
    spinup_day = 0.0,

    # Numerics
    damping_order = 4,
    damping_coef = 1.15741e-4,
    robert_coef = 0.04,
    implicit_coef = 0.5,

    # Restart
    # WARNING!!! Using a cold start would CLEANUP the restart directory!!!
    is_restart = false,
    restart_file = "",
    # is_restart = true,
    # restart_file = "/data92/garywu/undergrad_proposal/ctrl/restart/restart_t1728000000.jld2",
    saving_frequency = 86400 * 100,
    # saving_frequency = 0,    # disable saving restarts

    # Cold start (disabled if is_restart is true)
    initial_condition = :Moist_Spinup,

    # Physics
    moisture_processes = true,

    # IO
    output_path = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger = joinpath(output_path_base, "logger.log"),
    pressure_levels = [
        100000.0,
        92500.0,
        85000.0,
        70000.0,
        50000.0,
        30000.0,
        20000.0,
        10000.0,
        5000.0,
        1000.0,
    ],
    vars_to_output = [
        :u,
        :v,
        :w,
        :q,
        :t,
        :ps,
        :shflx,
        :lhflx,
        :precip,
        :bm_dt,
        :bm_dq,
        :bm_precip,
        :lrf_dt,
    ],
    output_interval = 21600,  # 6 hours

    do_plev_output = true,

    # Physics
    physics_params = physics_params,
)

# 4. Run Simulation
JGCM_Simulate(config)
