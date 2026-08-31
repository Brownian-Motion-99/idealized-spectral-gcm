using JGCM

# 1. Setup Paths
experiment_name = "Barotropic"
output_path_base = joinpath("exp", experiment_name)
mkpath(output_path_base)

# 2. Configure Model
config = Model_Config(
    name = experiment_name,
    model_type = :Barotropic,
    institution = "Group of Chaos and Predictability, Department of Atmospheric Sciences, National Taiwan University",

    # Resolution (T85 is standard for barotropic turbulence/instability tests)
    num_fourier = 85,
    nθ = 128,
    nd = 1,

    # Planet Settings
    radius = 6371.0e3,
    omega = 7.292e-5,
    grav = 9.80,
    day_to_sec = 86400,

    # Vertical Coordinate (None for Barotropic)
    # Using "even_sigma" with nd=1 is a safe fallback for the Driver logic
    vert_coord_option = "even_sigma",
    vert_difference_option = "none",
    vert_ref_level_option = "none",

    # Time Integration
    Δt = 600,       # Reduced Δt for stability at T85
    end_time = 86400 * 10, # 10-day simulation

    # Numerics
    damping_order = 4,
    damping_coef = 1.0e-4, # Hyperdiffusion to prevent spectral blocking
    robert_coef = 0.04,
    implicit_coef = 0.0,    # Barotropic is typically explicit

    # Restart logic
    is_restart        = false,
    restart_file      = "",
    saving_frequency = 86400,

    # Physics / Tracers
    moisture_processes = false,

    # IO Settings
    output_path = output_path_base,
    output_filename = joinpath(output_path_base, "barotropic_output.nc"),
    logger = joinpath(output_path_base, "logger.log"),

    # Variables must match symbols in Output_Mappings_Module for :Barotropic
    vars_to_output = [:u, :v, :vor, :ke],
    output_interval = 3600, # Hourly output
    do_plev_output  = false,

    # Initialization
    initial_condition = :Barotropic_Jet, # Standard Galewsky et al. or similar jet
    physics_params = Dict{String,Any}(),
)

# 3. Run Simulation
JGCM_Simulate(config)
