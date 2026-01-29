module Driver

using Printf
using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Dyn_Data_Module
using ..Vert_Coordinate_Module
using ..Time_Integrator_Module
using ..Output_Manager_Module
using ..Experiment_Configuration
using ..Initial_Conditions
using ..Semi_Implicit_Module
using ..Spectral_Dynamics_Module

using ..Barotropic_Dynamics_Module
using ..Shallow_Water_Dynamics_Module

export JGCM_Simulate

"""
    JGCM_Simulate(config::Model_Config)

The main entry point for running any JGCM simulation.
Handles allocation, initialization, time-stepping, and I/O based on the config.
"""
function JGCM_Simulate(config::Model_Config)
    
    msg_init       = "Initializing Experiment: $(config.name)"
    msg_model_info = "Model Type: $(config.model_type) | Resolution: T$(config.num_fourier)L$(config.nd)"

    @info msg_init
    @info msg_model_info
    
    open(config.logger, "w+") do log
        println(log, msg_init)
        println(log, msg_model_info)
    end

    # =========================================================================
    # 1. Core Allocation
    # =========================================================================
    
    # Mesh
    # Note: num_spherical usually = num_fourier + 1
    num_spherical = config.num_fourier + 1
    nλ = 2 * config.nθ
    
    mesh = Spectral_Spherical_Mesh(
        config.num_fourier, num_spherical, 
        nλ, config.nθ, config.nd, 
        config.radius
    )

    # 3D vs 2D Logic
    is_3d = (config.model_type == :PrimitiveEquation)
    
    # Vertical Coordinate (Only for 3D)
    vert_coord = nothing
    if is_3d
        # We assume standard sigma levels for now. 
        # Future: Add config.vert_coord_type to Model_Config for hybrid support.
        vert_coord = Vert_Coordinate(
            nλ, config.nθ, config.nd, 
            config.vert_coord_option, 
            config.vert_difference_option, 
            config.vert_ref_level_option
        )
    end

    # Atmo_Data
    # We parse physics flags from the config dictionary
    do_mass   = get(config.physics_params, "do_mass_correction",   true)
    do_energy = get(config.physics_params, "do_energy_correction", true)
    do_water  = get(config.physics_params, "do_water_correction",  false)
    use_virt  = get(config.physics_params, "use_virtual_temperature", false)

    atmo_data = Atmo_Data(
        config.name,
        nλ, config.nθ, config.nd,
        do_mass, do_energy, do_water, use_virt, config.L,
        mesh.sinθ;
        radius = config.radius, omega = config.omega, grav = config.grav,
        # Pass dictionary as kwargs to handle optional physics flags
        [Symbol(k) => v for (k,v) in config.physics_params]...
    )

    # Integrator
    # We map config.damping_coef etc. directly
    start_time = 0
    init_step  = true
    
    integrator = Filtered_Leapfrog(
        config.robert_coef, 
        config.damping_order, config.damping_coef, 
        mesh.laplacian_eig,
        config.implicit_coef,
        config.Δt, init_step, start_time, config.end_time
    )

    # Dyn_Data
    dyn_data = Dyn_Data(
        config.name, 
        config.num_fourier, num_spherical, 
        nλ, config.nθ, config.nd, 
        config.num_tracers
    )

    # Semi-Implicit Solver (Only for 3D)
    semi_implicit = nothing
    if is_3d && vert_coord !== nothing
        # Standard reference state
        ps_ref = 1.0e5
        t_ref  = fill(300.0, config.nd)
        semi_implicit = Semi_Implicit_Solver(vert_coord, atmo_data, integrator, ps_ref, t_ref, mesh.wave_numbers)
    end

    # =========================================================================
    # 2. Initialization
    # =========================================================================
    
    msg_init_cond = "Setting Initial Conditions..."
    @info msg_init_cond
    open(config.logger, "a") do log; println(log, msg_init_cond); end

    # Call the Bridge in Initial_Conditions.jl
    # Note: We pass full 'config' so it can access parameters (h0, perturbations)
    Initialize_Atmos_State!(mesh, atmo_data, dyn_data, vert_coord, config)
    
    # =========================================================================
    # 3. Output Management
    # =========================================================================
    
    # Map model type to Output_Manager symbol (:PrimitiveEquation, :Barotropic, etc.)
    # We reuse the logic from config.model_type but need to match Output_Manager expectations
    om_mode = if config.model_type == :PrimitiveEquation
        :PrimitiveEquation
    elseif config.model_type == :Shallow_Water
        :ShallowWater
    else
        :Barotropic
    end

    op_man = Output_Manager(
        mesh, vert_coord, atmo_data,
        start_time, config.end_time,
        config.vars_to_output;
        filename = config.output_filename,
        do_raw_output = config.do_raw_output,
        pressure_levels = config.pressure_levels,
        output_interval = config.output_interval,
        day_to_sec = config.day_to_sec,
        spinup_day = config.spinup_day,
        model_mode = om_mode
    )

    # Output Initial State
    Update_Output!(op_man, dyn_data, integrator.time)

    # =========================================================================
    # 4. Main Time Loop
    # =========================================================================
    
    NT = Int64(config.end_time / config.Δt)

    msg_start_loop = "Starting Time Loop: $NT steps"
    @info msg_start_loop
    open(config.logger, "a") do log; println(log, msg_start_loop); end
    
    # --- First Step (Euler / Init) ---
    Step_Dynamics!(
        config, mesh, atmo_data, dyn_data, 
        integrator, semi_implicit, vert_coord, config.physics_params
    )
    
    # If using Leapfrog, we need to correct the first step
    if isa(integrator, Filtered_Leapfrog)
        if isnothing(semi_implicit)
            Time_Integrator_Module.Update_Init_Step!(integrator)
        else
            Semi_Implicit_Module.Update_Init_Step!(semi_implicit)
        end
    end
    
    integrator.time += config.Δt
    Update_Output!(op_man, dyn_data, integrator.time)

    # --- Main Loop ---
    for i = 2:NT
        Step_Dynamics!(
            config, mesh, atmo_data, dyn_data, 
            integrator, semi_implicit, vert_coord, config.physics_params
        )
        
        integrator.time += config.Δt
        Update_Output!(op_man, dyn_data, integrator.time)

        # Simple Progress Log
        if i % (config.day_to_sec / config.Δt / 4) == 0
            
            day = integrator.time / config.day_to_sec

            msg_step_and_day = @sprintf("=== Step %d | Day %.2f ===", i, day)
            @info msg_step_and_day
            open(config.logger, "a") do log; println(log, msg_step_and_day); end
            
            # Define variables to monitor
            # Using Views for tracers to avoid allocation
            diag_vars = [
                (:U,      dyn_data.grid_u_c),
                (:V,      dyn_data.grid_v_c),
                (:T,      dyn_data.grid_t_c),
                (:W,      dyn_data.grid_w_full),
                (:P_full, dyn_data.grid_p_full),
                (:T_eq,   dyn_data.grid_t_eq),
                (:Q,      @view dyn_data.grid_tracers_c[:, :, :, 1])
            ]

            for (name, field) in diag_vars
                # Find index of maximum absolute magnitude
                max_val, idx = findmax(abs, field)
                
                # Retrieve the actual signed value at that location
                actual_val = field[idx]
                
                # Extract coordinates
                if ndims(field) == 3
                    
                    i, j, k = idx[1], idx[2], idx[3]
                    msg_var_diagnostic = @sprintf(
                        "%-6s: %10.4f  at (λ=%03d, θ=%03d, k=%02d)", 
                        string(name), actual_val, i, j, k
                    )
                    @info msg_var_diagnostic
                    open(config.logger, "a") do log; println(log, msg_var_diagnostic); end

                elseif ndims(field) == 2
                    
                    i, j = idx[1], idx[2]
                    msg_var_diagnostic = @sprintf(
                        "%-6s: %10.4f  at (λ=%03d, θ=%03d)", 
                        string(name), actual_val, i, j
                    )
                    @info msg_var_diagnostic
                    open(config.logger, "a") do log; println(log, msg_var_diagnostic); end

                end
            end
            
            println("-"^40) # Separator line
            open(config.logger, "a") do log; println(log, "-"^40); end

        end
    end
    
    Finalize_Output!(op_man)
    msg_end = "Simulation Complete."
    @info msg_end
    open(config.logger, "a") do log; println(log, msg_end); end

end

# ==============================================================================
# Helper: Dynamics Dispatcher
# ==============================================================================

function Step_Dynamics!(config, mesh, atmo_data, dyn_data, integrator, semi_implicit, vert_coord, physics_params)

    model_type = config.model_type
    if model_type == :Barotropic
        Barotropic_Dynamics!(mesh, atmo_data, dyn_data, integrator)
        
    elseif model_type == :Shallow_Water
        # Extract SW params safely
        kappa_m = get(physics_params, "kappa_m", 1.0/(20.0*86400))
        kappa_t = get(physics_params, "kappa_t", 1.0/(10.0*86400))
        h_eq    = dyn_data.grid_geopots # We stored h_eq here during init
        h_0     = get(physics_params, "h_0", 3.0e4)
        
        Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
        Shallow_Water_Dynamics!(mesh, atmo_data, h_0, dyn_data, integrator)
        
    elseif model_type == :PrimitiveEquation
        # Primitive Equation Dynamics
        # Note: Atmosphere_Update! handles physics internals based on Atmo_Data flags
        Atmosphere_Update!(config, mesh, vert_coord, atmo_data, dyn_data, semi_implicit, physics_params)
        
    else
        error("Unknown Model Type: $model_type")
    end
end

end