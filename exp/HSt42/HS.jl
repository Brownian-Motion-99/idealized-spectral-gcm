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


# function Atmos_Spectral_Dynamics_Main(physcis_params::Dict{String, Float64}, end_day::Int64 = 5, spinup_day::Int64 = 0, L::Float64 = L)
#     # the decay of a sinusoidal disturbance to a zonally symmetric flow 
#     # that resembles that found in the upper troposphere in Northern winter.
#     name = "Spectral_Dynamics"
#     num_fourier, nθ, nd = 42, 64, 20
#     #num_fourier, nθ, nd = 21, 32, 20

#     num_spherical = num_fourier + 1
#     nλ = 2nθ
    
#     radius = 6371000.0
#     omega = 7.292e-5
#     sea_level_ps_ref = 1.0e5
#     init_t = 264.0

    
#     # Initialize mesh
#     mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nθ, nλ, nd, radius)
#     θc, λc = mesh.θc,  mesh.λc
#     cosθ, sinθ = mesh.cosθ, mesh.sinθ
    
#     vert_coord = Vert_Coordinate(nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts", sea_level_ps_ref)
#     # Initialize atmo_data
#     do_mass_correction   = true
#     do_energy_correction = true
#     do_water_correction  = true
    
#     use_virtual_temperature = true
#     atmo_data = Atmo_Data(name, nλ, nθ, nd, do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature, sinθ, radius,  omega)
    
#     # Initialize integrator
#     damping_order = 8
#     damping_coef = 1.15741e-4
#     robert_coef  = 0.04 
    
#     implicit_coef = 0.5
#     day_to_sec = 86400
#     start_time = 0
#     end_time = end_day*day_to_sec  
#     Δt = 600
#     ### CJY
#     if warm_start_file_name != "None"
#         init_step = false # => In leapfrog would NOT do damping at initital time (should use in warm start case)
#     else
#         init_step = true # => In leapfrog would do damping at initital time
#     end
    
#     integrator = Filtered_Leapfrog(robert_coef, 
#     damping_order, damping_coef, mesh.laplacian_eig,
#     implicit_coef, Δt, init_step, start_time, end_time)
    
#     ps_ref = sea_level_ps_ref
#     t_ref = fill(300.0, nd)
#     wave_numbers = mesh.wave_numbers
#     semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data,
#     integrator, ps_ref, t_ref, wave_numbers)


#     # Data Visualization
#     op_man = Output_Manager(mesh, vert_coord, start_time, end_time, spinup_day)
        
    
#     # Initialize data
#     # By CJY edit for passive tracer
#     num_grid_tracters = 1
#     num_spe_tracters  = 1
#     dyn_data = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd,num_grid_tracters ,num_spe_tracters) ### origin = Dyn_Data(name, num_fourier, num_spherical, nλ, nθ, nd)

#     NT = Int64(end_time / Δt)
    
#     Get_Topography!(dyn_data.grid_geopots, warm_start_file_name, initial_day)
    
#     Spectral_Initialize_Fields!(mesh, atmo_data, vert_coord, sea_level_ps_ref, init_t, dyn_data.grid_geopots, dyn_data.T_ref, dyn_data, Δt, warm_start_file_name, initial_day)
    
    
#     Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data, physcis_params, L, dyn_data.T_ref)
#     Update_Init_Step!(semi_implicit)
#     integrator.time += Δt
#     Update_Output!(op_man, dyn_data, integrator.time)
    
    
#     for i = 2:NT

#         Atmosphere_Update!(mesh, atmo_data, vert_coord, semi_implicit, dyn_data, physcis_params, L, dyn_data.T_ref)
#         integrator.time += Δt
#         #@info integrator.time

#         Update_Output!(op_man, dyn_data, integrator.time)

#         # if (integrator.time%day_to_sec == 0)
#         #     # dyn_data.grid_tracers_c[dyn_data.grid_tracers_c .< 0] .= 0
#         #     @info "Day: ", div(integrator.time,day_to_sec), " Max |U|,|V|,|P|,|T|,|qv|: ", maximum(abs.(dyn_data.grid_u_c)), maximum(abs.(dyn_data.grid_v_c)), maximum(dyn_data.grid_p_full), maximum(dyn_data.grid_t_c), maximum(dyn_data.grid_tracers_c)
#         #     @info "Day: ", div(integrator.time,day_to_sec), " Min |U|,|V|,|P|,|T|,|qv|: ", minimum(abs.(dyn_data.grid_u_c)), minimum(abs.(dyn_data.grid_v_c)), minimum(dyn_data.grid_p_full), minimum(dyn_data.grid_t_c), minimum(dyn_data.grid_tracers_c)
#         # end
#     end

#     return op_man
    
# end

