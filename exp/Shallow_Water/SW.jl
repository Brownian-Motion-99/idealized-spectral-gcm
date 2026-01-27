using JGCM

using JGCM

# 1. Define Physics Parameters
# These control the forcing and damping specific to the SW test case
physics_params = Dict{String, Any}(
    "h_0"        => 3.0e4,    # Reference depth (m)
    "h_amp"      => 2.0e4,    # Amplitude of the forcing bump (m)
    "h_lon"      => 90.0,     # Longitude of forcing center (deg)
    "h_lat"      => 25.0,     # Latitude of forcing center (deg)
    "h_width"    => 15.0,     # Width of forcing Gaussian (deg)
    "alpha"      => 0.0,      # Rotation angle (if needed for rotated pole)
    "h_itcz"     => 1.0e5,    # ITCZ Amplitude
    "itcz_width" => 4.0,      # ITCZ Width (deg)
    
    # Damping timescales (Newtonian cooling & Rayleigh friction)
    "kappa_m" => 1.0/(20.0*86400.0), # Momentum damping (20 days)
    "kappa_t" => 1.0/(10.0*86400.0)  # Thermal relaxation (10 days)
)

# 2. Configure the Model
config = Model_Config(
    name = "Shallow_Water_Test",
    model_type = :Shallow_Water,
    
    # Resolution
    num_fourier = 85,
    nθ = 128, 
    nd = 1,
    
    # Planet
    radius = 6371.0e3, 
    omega = 7.292e-5,
    grav = 9.80,

    vert_coord_option = nothing, vert_difference_option = nothing, vert_ref_level_option = nothing,
    
    # Time Integration
    day_to_sec = 86400,
    Δt = 1200, 
    end_time = 86400 * 20,  # Run for 20 days
    
    # Numerics
    damping_order = 4, 
    damping_coef = 1.0e-4, 
    robert_coef = 0.04, 
    implicit_coef = 0.5, # Important for gravity waves!
    
    # IO
    output_filename = "exp/Shallow_Water/output.nc",
    logger = "exp/Shallow_Water/logger.log",
    vars_to_output = [:u, :v, :h, :vor, :div, :pv], # :pv = Potential Vorticity
    output_interval = 1200,
    
    # Initialization & Physics
    num_grid_tracters = 1,
    num_spe_tracters  = 1,
    initial_condition = :Shallow_Water_Test,
    physics_params = physics_params
)

# 3. Run
JGCM_Simulate(config)