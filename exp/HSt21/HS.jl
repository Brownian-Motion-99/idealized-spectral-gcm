using JGCM

# 1. Physics Configuration
physics_params = Dict{String, Any}(
    
    # Corrections    
    "do_mass_correction"   => true,
    "do_energy_correction" => true,
    "do_water_correction"  => true,

    # Grid scale condensation
    "do_Lscale_Cond" => false,
    "L"              => 0.2,

    # PBL fluxes
    "do_Sensible_Heating"    => false,
    "C_H"                    => 0.0044,
    "do_Surface_Evaporation" => false,
    "C_E"                    => 0.0044,
    "do_Implicit_PBL_Scheme" => false,
    "C_D"                    => 0.0044,
    
    # "PBL_Top_Mode"  => :ModelLevel,
    # "PBL_Top_Value" => 4,
    "PBL_Top_Mode"  => :PressureLevel,
    "PBL_Top_Value" => 85000.0,

    # Held-Suarez
    "do_HS_Forcing" => true,
    "σ_b"           => 0.7,
    "k_a"           => 1.0/(40.0),
    "k_s"           => 1.0/(4.0),
    "k_f"           => 1.0/(1.0),
    "ΔT_y"          => 60.0, 
    "Δθ_z"          => 10.0

    # TODO: cumulus parameterization
    # TODO: radiation parameterization
)

# 2. Define Output Paths *Before* Configuration
experiment_name  = "HSt21"
output_path_base = joinpath("exp", experiment_name) 
mkpath(output_path_base)

# 3. Model Configuration
config = Model_Config(

    name = "HS_Moist_T21",
    model_type = :PrimitiveEquation,
    
    # Resolution
    num_fourier = 21, 
    nθ          = 32, 
    nd          = 10, 
    
    # Vertical Coordinate
    vert_coord_option      = "even_sigma", 
    vert_difference_option = "simmons_and_burridge", 
    vert_ref_level_option  = "second_centered_wts",
    
    # Planet Settings
    radius     = 6371.0e3, 
    omega      = 7.292e-5, 
    grav       = 9.80,
    day_to_sec = 86400,
    
    # Time Integration
    Δt         = 600,
    end_time   = 86400 * 4,
    spinup_day = 0.0,
    
    # Numerics
    damping_order = 4, 
    damping_coef  = 1.15741e-4, 
    robert_coef   = 0.04, 
    implicit_coef = 0.5,
    
    # Restart
    # WARNING!!! Using a cold start would CLEANUP the restart directory!!!
    is_restart        = false,
    restart_file      = "",
    # is_restart        = true,
    # restart_file      = joinpath(output_path_base, "restart", "restart_t86400.jld2"),
    # restart_frequency = 86400 * 2,
    restart_frequency = 0,    # disable saving restarts

    # Cold start (disabled if is_restart is true)
    initial_condition = :Moist_Spinup,

    # Physics
    moisture_processes = true,
    num_tracers = 1,
    
    # IO
    output_path     = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger          = joinpath(output_path_base, "logger.log"),

    do_raw_output   = true,
    pressure_levels = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    vars_to_output  = [:u, :v, :w, :q, :t, :ps, :shflx, :lhflx, :precip],
    output_interval = 600,
    
    # Physics
    physics_params = physics_params
)

# 4. Run Simulation
JGCM_Simulate(config)