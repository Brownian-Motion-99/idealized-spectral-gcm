using JGCM

config = Model_Config(
    name = "Barotropic",
    model_type = :Barotropic,

    num_fourier = 85, nθ = 128, nd = 1,
    radius = 6371.0e3, omega = 7.292e-5, grav = 9.80,
    
    vert_coord_option = nothing, vert_difference_option = nothing, vert_ref_level_option = nothing,
    
    Δt = 1800, end_time = 691200, day_to_sec = 86400,
    
    damping_order = 4, damping_coef = 1.e-4, 
    robert_coef = 0.04, implicit_coef = 0.0,

    num_tracers = 0, 
    
    initial_condition = :Barotropic_Jet,
    
    output_filename = "exp/Barotropic/output_test.nc",
    logger = "exp/Barotropic/logger.log",
    vars_to_output = [:u, :v, :vor],
    output_interval = 1800,

    physics_params = Dict{String, Any}()
)

JGCM_Simulate(config)