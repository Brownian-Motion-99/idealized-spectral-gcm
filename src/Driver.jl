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
using ..Restart_Manager_Module

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
    do_water  = get(config.physics_params, "do_water_correction",  true)
    use_virt  = get(config.physics_params, "use_virtual_temperature", true)

    atmo_data = Atmo_Data(
        config.name,
        nλ, config.nθ, config.nd,
        do_mass, do_energy, do_water, use_virt,
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
    
    # Setup Restart Manager
    saving_freq = config.saving_frequency
    restart_dir = joinpath(config.output_path, "restart")
    restart_mgr = Restart_Manager(restart_dir, saving_freq)

    # State Variables for Integrator
    start_time = 0
    init_step  = true

    if config.is_restart
        # --- PATH A: WARM START ---
        msg_restart = "Warm Start Detected. Loading: $(config.restart_file)"
        @info msg_restart
        open(config.logger, "a") do log; println(log, msg_restart); end

        if !isfile(config.restart_file)
            error("Warm start requested but file missing: $(config.restart_file)")
        end

        # Load data AND get the time we left off
        saved_time = Load_Restart_File!(dyn_data, config.restart_file)

        start_time = saved_time
        init_step  = false  # We are resuming, so we don't need the Euler start

        # Sync integrator to the correct restart time (integrator was constructed
        # before the restart file was loaded, so its time fields are still 0)
        integrator.time       = start_time
        integrator.start_time = start_time
        integrator.init_step  = false
        
    else
        # --- PATH B: COLD START ---
        msg_init_cond = "Cold Start: Setting Analytical Initial Conditions..."
        @info msg_init_cond
        open(config.logger, "a") do log; println(log, msg_init_cond); end

        if isdir(restart_dir) && !isempty(readdir(restart_dir))
            msg_clean_restart = "COLD START DETECTED: Cleaning up old restart files in $restart_dir"
            @warn msg_clean_restart
            open(config.logger, "a") do log; println(log, msg_clean_restart); end
            for f in readdir(restart_dir, join=true)
                if endswith(f, ".jld2")
                    rm(f)
                end
            end
        end

        # Standard initialization
        Initialize_Atmos_State!(mesh, atmo_data, dyn_data, vert_coord, config)
        
        start_time = 0
        init_step  = true
    end
    
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

    # Sanity check: do_plev_output requires pressure_levels to be set
    if config.do_plev_output && isempty(config.pressure_levels)
        error("do_plev_output = true but pressure_levels is empty. Please specify the target pressure levels.")
    end

    op_man = Output_Manager(
        mesh, vert_coord, atmo_data,
        start_time, config.end_time,
        config.vars_to_output;
        filename        = config.output_filename,
        do_plev_output  = config.do_plev_output,
        pressure_levels = config.pressure_levels,
        output_interval = config.output_interval,
        day_to_sec      = config.day_to_sec,
        spinup_day      = config.spinup_day,
        model_mode      = om_mode,
        institute       = config.institution,
        experiment_id   = config.name
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

    msg_cleanup_old_restarts = "Only the last 5 restart files are kept to save space."
    @warn msg_cleanup_old_restarts
    open(config.logger, "a") do log; println(log, msg_cleanup_old_restarts); end
    
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

        # Checkpoint + NC chunk rotation (coordinated at saving_frequency)
        if restart_mgr.saving_frequency > 0 && integrator.time > 0 && (integrator.time % restart_mgr.saving_frequency == 0)
            Write_Restart_File(restart_mgr, dyn_data, Int64(integrator.time))

            # Don't rotate to a fresh NC chunk on the final step: the chunk just
            # flushed above already holds the last interval's data, and a new
            # chunk opened here would never receive any (the loop ends next).
            if integrator.time < config.end_time
                Rotate_NC_Chunk!(op_man, Int64(integrator.time))
            end

            msg_ckpt = "Checkpoint saved at t=$(integrator.time)"
            @info msg_ckpt
            open(config.logger, "a") do log; println(log, msg_ckpt); end

            # Keep only the last 5 files to save space
            # Keep the starting file
            if config.is_restart
                Restart_Manager_Module.Cleanup_Old_Restarts(restart_mgr, 5, config.restart_file)
            else
                Restart_Manager_Module.Cleanup_Old_Restarts(restart_mgr, 5)
            end
        end

        # Simple Progress Log
        if i % (config.day_to_sec / config.Δt / 4) == 0
            status_diagnostics(i, config, dyn_data, integrator)
        end
    end
    
    Finalize_Output!(op_man)
    msg_end = "Simulation Complete."
    @info msg_end
    open(config.logger, "a") do log; println(log, msg_end); end

end



# ==============================================================================
# Helper: Status Diagnostics
# ==============================================================================
function status_diagnostics(step::Int64, config::Model_Config, dyn_data::Dyn_Data, integrator::Filtered_Leapfrog)

    day = integrator.time / config.day_to_sec

    msg_step_and_day = @sprintf("=== Step %d | Day %.2f ===", step, day)
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
        (:Q,      dyn_data.grid_q_c[:, :, :, 1])
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