using JGCM

# 1. Setup Paths
experiment_name = "Barotropic"
output_path_base = joinpath("exp", experiment_name)
mkpath(output_path_base)

# 2. Configure Model
config = Model_Config(
    name = experiment_name,
    model_type = :Barotropic,

    # Resolution
    num_fourier = 85, nθ = 128, nd = 1,
    
    # Planet
    radius = 6371.0e3, omega = 7.292e-5, grav = 9.80,
    
    # Vertical Coordinate (None for Barotropic)
    vert_coord_option = nothing, 
    vert_difference_option = nothing, 
    vert_ref_level_option = nothing,
    
    # Time Integration
    Δt = 1800, end_time = 691200, day_to_sec = 86400,
    
    # Numerics
    damping_order = 4, damping_coef = 1.e-4, robert_coef = 0.04, implicit_coef = 0.0,

    # Restart
    is_restart = false,
    restart_file = "",
    restart_frequency = 86400,

    # Physics / Tracers
    L = 0.0,
    moisture_processes = false,
    num_tracers = 0, # Usually 0 for pure barotropic
    
    # IO
    output_path     = output_path_base,
    output_filename = joinpath(output_path_base, "output_test.nc"),
    logger          = joinpath(output_path_base, "logger.log"),
    
    vars_to_output  = [:u, :v, :vor],
    output_interval = 1800,

    # Initialization
    initial_condition = :Barotropic_Jet,
    physics_params    = Dict{String, Any}()
)

JGCM_Simulate(config)