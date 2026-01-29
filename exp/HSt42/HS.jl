using JGCM

# 1. Physics Configuration
# These boolean flags control the complexity of the physics package.
# For a "Dry" Held-Suarez run, set do_moist_phys = false.
physics_params = Dict{String, Any}(
    # --- Dynamical Core Corrections ---
    "do_mass_correction"   => true,
    "do_energy_correction" => true,
    "do_water_correction"  => true,
    
    # --- Physical Parameterizations ---
    "do_large_scale_condensation" => true,
    "do_Sensible_heat_fluxes"     => true,
    "do_Surface_evaporation"      => true,
    "do_Implicit_PBL_Scheme"      => true,
    
    # --- Forcing Parameters (Held-Suarez) ---
    "σ_b"     => 0.7,
    "k_a"     => 1.0/(40.0),
    "k_s"     => 1.0/(4.0),
    "k_f"     => 1.0/(1.0),
    "ΔT_y"    => 60.0, 
    "Δθ_z"    => 10.0
)

# 2. Model Configuration
config = Model_Config(

    name = "HS_Moist_T42",
    model_type = :PrimitiveEquation, # 3D Mode
    
    # Resolution (T42L20)
    num_fourier = 42, 
    nθ          = 64, 
    nd          = 20, 
    
    # Vertical Coordinate
    vert_coord_option = "even_sigma", vert_difference_option = "simmons_and_burridge", vert_ref_level_option = "second_centered_wts",
    
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
    
    # L (latent heating parameter)
    L = 0.2,

    # Tracers (Important for Moist runs!)
    num_grid_tracters = 1,
    num_spe_tracters  = 1,
    
    # IO
    output_filename = "exp/HSt42/output.nc",
    logger = "exp/HSt42/logger.log",
    do_raw_output = true,
    pressure_levels = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    vars_to_output = [:u, :v, :w, :q, :t, :ps, :shflx, :lhflx],
    output_interval = 600,
    
    # Initialization
    initial_condition = :Moist_Spinup, # Calls Init_3D_Standard!
    physics_params = physics_params

)

# 3. Run Simulation
JGCM_Simulate(config)