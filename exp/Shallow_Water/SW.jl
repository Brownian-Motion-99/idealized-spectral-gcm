using JGCM

# 1. Define Physics Parameters
physics_params = Dict{String,Any}(
    "h_0" => 3.0e4,
    "h_amp" => 2.0e4,
    "h_lon" => 90.0,
    "h_lat" => 25.0,
    "h_width" => 15.0,
    "alpha" => 0.0,
    "h_itcz" => 1.0e5,
    "itcz_width" => 4.0,
    "kappa_m" => 1.0 / (20.0 * 86400.0),
    "kappa_t" => 1.0 / (10.0 * 86400.0),
)

# 2. Setup Paths
experiment_name = "Shallow_Water"
output_path_base = joinpath("exp", experiment_name)
mkpath(output_path_base)

# 3. Configure the Model
config = Model_Config(
    name = experiment_name,
    model_type = :Shallow_Water,
    institution = "Group of Chaos and Predictability, Department of Atmospheric Sciences, National Taiwan University",

    # Resolution
    num_fourier = 85,
    nθ = 128,
    nd = 1,

    # Planet
    radius = 6371.0e3,
    omega = 7.292e-5,
    grav = 9.80,

    # Vertical Coordinate (None for SW)
    vert_coord_option = nothing,
    vert_difference_option = nothing,
    vert_ref_level_option = nothing,

    # Time Integration
    day_to_sec = 86400,
    Δt = 1200,
    end_time = 86400 * 20,

    # Numerics
    damping_order = 4,
    damping_coef = 1.0e-4,
    robert_coef = 0.04,
    implicit_coef = 0.5,

    # Restart (Default Off)
    is_restart = false,
    restart_file = "",
    saving_frequency = 86400,

    # Physics / Tracers
    moisture_processes = false,
    num_tracers = 1, # Kept as 1 for potential vorticity tracer if needed

    # IO
    output_path = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger = joinpath(output_path_base, "logger.log"),
    vars_to_output = [:u, :v, :h, :vor, :div, :pv],
    output_interval = 1200,

    # Initialization
    initial_condition = :Shallow_Water_Test,
    physics_params = physics_params,
)

# 4. Run
JGCM_Simulate(config)
