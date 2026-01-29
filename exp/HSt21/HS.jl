using JGCM

# 1. Physics Configuration
physics_params = Dict{String, Any}(
    "do_mass_correction"   => true,
    "do_energy_correction" => true,
    "do_water_correction"  => true,
    "do_Lscale_Cond"          => true,
    "do_Sensible_Heating"     => true,
    "do_Surface_Evaporation"  => true,
    "do_Implicit_PBL_Scheme"  => true,
    
    "PBL_Top_Mode"  => :ModelLevel,
    "PBL_Top_Value" => 2,
    # "PBL_Top_Mode"  => :PressureLevel,
    # "PBL_Top_Value" => 85000.0,

    "σ_b"     => 0.7,
    "k_a"     => 1.0/(40.0),
    "k_s"     => 1.0/(4.0),
    "k_f"     => 1.0/(1.0),
    "ΔT_y"    => 60.0, 
    "Δθ_z"    => 10.0
)

# 2. Setup Paths
experiment_name = "HSt21" # Changed folder name to match T21 context
output_path_base = joinpath("exp", experiment_name)
mkpath(output_path_base)

# 3. Model Configuration
config = Model_Config(

    name = "HS_Moist_T21",
    model_type = :PrimitiveEquation,
    
    # Resolution (T21L10)
    num_fourier = 21, nθ = 32, nd = 10, 
    
    # Vertical Coordinate
    vert_coord_option = "even_sigma", 
    vert_difference_option = "simmons_and_burridge", 
    vert_ref_level_option = "second_centered_wts",
    
    # Planet Settings
    radius = 6371.0e3, omega = 7.292e-5, grav = 9.80, day_to_sec = 86400,
    
    # Time Integration
    Δt = 600, end_time = 86400 * 4, spinup_day = 0.0,
    
    # Numerics
    damping_order = 4, damping_coef = 1.15741e-4, robert_coef = 0.04, implicit_coef = 0.5,
    
    # Restart
    is_restart = false,
    restart_file = "",
    restart_frequency = 86400,

    # Physics
    L = 0.2,
    moisture_processes = true,
    num_tracers = 1, # Renamed from num_grid_tracters/num_spe_tracters
    
    # IO
    output_path     = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger          = joinpath(output_path_base, "logger.log"),
    
    do_raw_output   = true,
    pressure_levels = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    vars_to_output  = [:u, :v, :w, :q, :t, :ps, :shflx, :lhflx],
    output_interval = 600,
    
    # Initialization
    initial_condition = :Moist_Spinup,
    physics_params    = physics_params
)

# 4. Run Simulation
JGCM_Simulate(config)