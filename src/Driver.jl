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
using ..Atmos_Param_Module

using ..Spectral_Dynamics_Module
using ..Barotropic_Dynamics_Module
using ..Shallow_Water_Dynamics_Module

export JGCM_Simulate

"""
    remaining_time_steps(start_time, end_time, Δt)

Return the number of timesteps needed to advance from `start_time` to
`end_time`. The interval must be nonnegative and exactly divisible by `Δt`.
"""
function remaining_time_steps(start_time::Integer, end_time::Integer, Δt::Integer)
    Δt > 0 || throw(ArgumentError("Δt must be positive, got $Δt"))
    end_time >= start_time || throw(
        ArgumentError(
            "end_time ($end_time) must not be earlier than start_time ($start_time)",
        ),
    )

    remaining_time = end_time - start_time
    remaining_time % Δt == 0 || throw(
        ArgumentError(
            "remaining simulation time ($remaining_time) must be evenly divisible by Δt ($Δt)",
        ),
    )

    return Int64(div(remaining_time, Δt))
end

function progress_metrics(current_time, start_time, end_time, Δt, elapsed_seconds)
    current_time <= end_time || throw(
        ArgumentError("current_time ($current_time) exceeds end_time ($end_time)"),
    )
    completed_steps = remaining_time_steps(start_time, current_time, Δt)
    total_steps = remaining_time_steps(start_time, end_time, Δt)
    segment_progress = total_steps == 0 ? 1.0 : completed_steps / total_steps
    overall_progress = end_time == 0 ? 1.0 : current_time / end_time
    eta_seconds = completed_steps == 0 ? nothing :
                  elapsed_seconds * (total_steps - completed_steps) / completed_steps

    return (; completed_steps, total_steps, segment_progress, overall_progress, eta_seconds)
end

function run_metrics(start_time, current_time, completed_steps, elapsed_seconds, day_to_sec)
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    simulated_days = (current_time - start_time) / day_to_sec
    seconds_per_step = completed_steps == 0 ? nothing : elapsed_seconds / completed_steps
    simulated_days_per_wall_day = completed_steps == 0 || elapsed_seconds <= 0 ? nothing :
                                  simulated_days * 86400 / elapsed_seconds

    return (; simulated_days, seconds_per_step, simulated_days_per_wall_day)
end

"""
    JGCM_Simulate(config::Model_Config)

The main entry point for running any JGCM simulation.
Handles allocation, initialization, time-stepping, and I/O based on the config.
"""
function JGCM_Simulate(config::Model_Config)
    simulation_start_ns = time_ns()
    
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
        config.num_fourier,
        num_spherical,
        nλ,
        config.nθ,
        config.nd,
        config.radius,
    )

    # 3D vs 2D Logic
    is_3d = (config.model_type == :PrimitiveEquation)

    # Vertical Coordinate (Only for 3D)
    vert_coord = nothing
    if is_3d
        vert_coord = Vert_Coordinate(
            nλ,
            config.nθ,
            config.nd,
            config.vert_coord_option,
            config.vert_difference_option,
            config.vert_ref_level_option,
        )
    end

    # Atmo_Data
    # We parse physics flags from the config dictionary
    do_mass = get(config.physics_params, "do_mass_correction", true)
    do_energy = get(config.physics_params, "do_energy_correction", true)
    do_water = get(config.physics_params, "do_water_correction", true)
    use_virt = get(config.physics_params, "use_virtual_temperature", true)

    atmo_data = Atmo_Data(
        config.name,
        nλ,
        config.nθ,
        config.nd,
        do_mass,
        do_energy,
        do_water,
        use_virt,
        mesh.sinθ;
        radius = config.radius,
        omega = config.omega,
        grav = config.grav,
        # Pass dictionary as kwargs to handle optional physics flags
        [Symbol(k) => v for (k, v) in config.physics_params]...,
    )

    # Integrator
    # We map config.damping_coef etc. directly
    start_time = 0
    init_step = true

    integrator = Filtered_Leapfrog(
        config.robert_coef,
        config.damping_order,
        config.damping_coef,
        mesh.laplacian_eig,
        config.implicit_coef,
        config.Δt,
        init_step,
        start_time,
        config.end_time,
    )

    # Dyn_Data
    dyn_data = Dyn_Data(
        config.name,
        config.num_fourier,
        num_spherical,
        nλ,
        config.nθ,
        config.nd,
        config.num_tracers,
    )

    # Load static LRF data once, before the time loop. The retained state is
    # stored in the existing physics configuration dictionary for now.
    if get(config.physics_params, "do_LRF", false)
        config.moisture_processes || error("LRF requires moisture_processes = true")
        lrf_file = get(config.physics_params, "LRF_file", nothing)
        lrf_file isa AbstractString || error("LRF requires physics_params[\"LRF_file\"]")
        config.physics_params["LRF_state"] =
            Load_LRF_State(lrf_file, nλ, config.nθ, config.nd)
    end

    # Construct Betts-Miller configuration and reusable column work arrays once.
    if get(config.physics_params, "do_Betts_Miller", false)
        config.moisture_processes ||
            error("Betts-Miller requires moisture_processes = true")
        bm_tau = Float64(get(config.physics_params, "bm_tau", 7200.0))
        bm_relative_humidity =
            Float64(get(config.physics_params, "bm_relative_humidity", 0.8))
        2 * config.Δt <= bm_tau || throw(
            ArgumentError(
                "Betts-Miller requires 2 * Δt <= bm_tau for leapfrog integration; " *
                "got Δt=$(config.Δt) s and bm_tau=$bm_tau s",
            ),
        )
        config.physics_params["BM_state"] = Betts_Miller_State(
            config.nd;
            tau = bm_tau,
            relative_humidity = bm_relative_humidity,
        )
    end

    # Semi-Implicit Solver (Only for 3D)
    semi_implicit = nothing
    if is_3d && vert_coord !== nothing
        # Standard reference state
        ps_ref = 1.0e5
        t_ref = fill(300.0, config.nd)
        semi_implicit = Semi_Implicit_Solver(
            vert_coord,
            atmo_data,
            integrator,
            ps_ref,
            t_ref,
            mesh.wave_numbers,
        )
    end

    # =========================================================================
    # 2. Initialization
    # =========================================================================

    # Setup Restart Manager
    restart_freq = config.restart_frequency
    restart_dir = joinpath(config.output_path, "restart")
    restart_mgr = Restart_Manager(restart_dir, restart_freq)

    # State Variables for Integrator
    start_time = 0
    init_step = true

    if config.is_restart
        # --- PATH A: WARM START ---
        msg_restart = "Warm Start Detected. Loading: $(config.restart_file)"
        @info msg_restart
        open(config.logger, "a") do log
            println(log, msg_restart)
        end

        if !isfile(config.restart_file)
            error("Warm start requested but file missing: $(config.restart_file)")
        end

        # Load data AND get the time we left off
        saved_time = Load_Restart_File!(dyn_data, config.restart_file)

        start_time = saved_time
        init_step = false  # We are resuming, so we don't need the Euler start

        # Sync integrator to the correct restart time (integrator was constructed
        # before the restart file was loaded, so its time fields are still 0)
        integrator.time = start_time
        integrator.start_time = start_time
        integrator.init_step  = false

        # The semi-implicit solver was constructed with the cold-start timestep.
        # Rebuild its wave matrix for the leapfrog timestep used after a restart.
        if !isnothing(semi_implicit)
            Semi_Implicit_Module.Update_Init_Step!(semi_implicit)
        end
        
    else
        # --- PATH B: COLD START ---
        msg_init_cond = "Cold Start: Setting Analytical Initial Conditions..."
        @info msg_init_cond
        open(config.logger, "a") do log
            println(log, msg_init_cond)
        end

        if isdir(restart_dir) && !isempty(readdir(restart_dir))
            msg_clean_restart = "COLD START DETECTED: Cleaning up old restart files in $restart_dir"
            @warn msg_clean_restart
            open(config.logger, "a") do log
                println(log, msg_clean_restart)
            end
            for f in readdir(restart_dir, join = true)
                if endswith(f, ".jld2")
                    rm(f)
                end
            end
        end

        # Standard initialization
        Initialize_Atmos_State!(mesh, atmo_data, dyn_data, vert_coord, config)

        start_time = 0
        init_step = true
    end

    NT = remaining_time_steps(start_time, config.end_time, config.Δt)
    
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
        mesh,
        vert_coord,
        atmo_data,
        start_time,
        config.end_time,
        config.vars_to_output;
        filename = config.output_filename,
        do_raw_output = config.do_raw_output,
        pressure_levels = config.pressure_levels,
        output_interval = config.output_interval,
        day_to_sec = config.day_to_sec,
        spinup_day = config.spinup_day,
        model_mode = om_mode,
        institute = config.institution,
        experiment_id = config.name,
    )

    # Output Initial State
    Update_Output!(op_man, dyn_data, integrator.time)

    # =========================================================================
    # 4. Main Time Loop
    # =========================================================================
    
    msg_start_loop = "Starting Time Loop: $NT remaining steps"
    @info msg_start_loop
    open(config.logger, "a") do log
        println(log, msg_start_loop)
    end

    msg_cleanup_old_restarts = "Only the last 5 restart files are kept to save space."
    @warn msg_cleanup_old_restarts
    open(config.logger, "a") do log
        println(log, msg_cleanup_old_restarts)
    end

    loop_start_ns = time_ns()

    if NT > 0
        # --- First Step (Euler / Init for a cold start) ---
        Step_Dynamics!(
            config, mesh, atmo_data, dyn_data,
            integrator, semi_implicit, vert_coord, config.physics_params
        )

        # A warm restart has init_step = false and continues with leapfrog directly.
        if integrator.init_step
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

            # Restart
            if restart_mgr.restart_frequency > 0 && integrator.time > 0 && (integrator.time % restart_mgr.restart_frequency == 0)
                Write_Restart_File(restart_mgr, dyn_data, Int64(integrator.time))

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

            # Simple progress indicator
            if i % (config.day_to_sec / config.Δt / 4) == 0
                elapsed_seconds = (time_ns() - loop_start_ns) / 1.0e9
                status_diagnostics(elapsed_seconds, config, dyn_data, integrator)
            end
        end
    end

    Finalize_Output!(op_man)
    integration_seconds = (time_ns() - loop_start_ns) / 1.0e9
    initialization_seconds = (loop_start_ns - simulation_start_ns) / 1.0e9
    total_seconds = (time_ns() - simulation_start_ns) / 1.0e9
    completed_steps = remaining_time_steps(start_time, integrator.time, config.Δt)
    metrics = run_metrics(
        start_time,
        integrator.time,
        completed_steps,
        integration_seconds,
        config.day_to_sec,
    )

    msg_end = "Simulation Complete."
    msg_metrics = if completed_steps == 0
        @sprintf(
            "Run Metrics: Model Day %.2f -> %.2f | Steps: 0 | Initialization: %s | Integration: %s | Total: %s | Seconds/Step: N/A | Simulated Days/Wall Day: N/A",
            start_time / config.day_to_sec,
            integrator.time / config.day_to_sec,
            format_duration(initialization_seconds),
            format_duration(integration_seconds),
            format_duration(total_seconds),
        )
    else
        @sprintf(
            "Run Metrics: Model Day %.2f -> %.2f | Steps: %d | Initialization: %s | Integration: %s | Total: %s | %.3f s/step | %.2f simulated days/wall day",
            start_time / config.day_to_sec,
            integrator.time / config.day_to_sec,
            completed_steps,
            format_duration(initialization_seconds),
            format_duration(integration_seconds),
            format_duration(total_seconds),
            metrics.seconds_per_step,
            metrics.simulated_days_per_wall_day,
        )
    end
    @info msg_end
    @info msg_metrics
    open(config.logger, "a") do log
        println(log, msg_end)
        println(log, msg_metrics)
    end

end



# ==============================================================================
# Helper: Status Diagnostics
# ==============================================================================
function status_diagnostics(
    elapsed_seconds::Float64,
    config::Model_Config,
    dyn_data::Dyn_Data,
    integrator::Filtered_Leapfrog,
)

    day = integrator.time / config.day_to_sec
    metrics = progress_metrics(
        integrator.time,
        integrator.start_time,
        config.end_time,
        config.Δt,
        elapsed_seconds,
    )

    msg_step_and_day = @sprintf(
        "=== Segment Step %d/%d | Day %.2f ===",
        metrics.completed_steps,
        metrics.total_steps,
        day,
    )
    @info msg_step_and_day
    open(config.logger, "a") do log
        println(log, msg_step_and_day)
    end

    msg_progress = @sprintf(
        "Segment: %.1f%% | Overall: %.1f%% | Elapsed: %s | ETA: %s",
        100 * metrics.segment_progress,
        100 * metrics.overall_progress,
        format_duration(elapsed_seconds),
        isnothing(metrics.eta_seconds) ? "N/A" : format_duration(metrics.eta_seconds),
    )
    @info msg_progress
    open(config.logger, "a") do log
        println(log, msg_progress)
    end

    # Define variables to monitor
    # Using Views for tracers to avoid allocation
    diag_vars = [
        (:U, dyn_data.grid_u_c),
        (:V, dyn_data.grid_v_c),
        (:T, dyn_data.grid_t_c),
        (:W, dyn_data.grid_w_full),
        (:P_full, dyn_data.grid_p_full),
        (:T_eq, dyn_data.grid_t_eq),
        (:Q, dyn_data.grid_q_c[:, :, :, 1]),
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
                string(name),
                actual_val,
                i,
                j,
                k
            )
            @info msg_var_diagnostic
            open(config.logger, "a") do log
                println(log, msg_var_diagnostic)
            end

        elseif ndims(field) == 2

            i, j = idx[1], idx[2]
            msg_var_diagnostic = @sprintf(
                "%-6s: %10.4f  at (λ=%03d, θ=%03d)",
                string(name),
                actual_val,
                i,
                j
            )
            @info msg_var_diagnostic
            open(config.logger, "a") do log
                println(log, msg_var_diagnostic)
            end

        end
    end

    println("-"^40) # Separator line
    open(config.logger, "a") do log
        println(log, "-"^40)
    end

end


"""
Format a duration in a readable form for progress log.
"""
function format_duration(seconds::Real)
    total_seconds = max(0, round(Int, seconds))
    days, remainder = divrem(total_seconds, 86400)
    hours, remainder = divrem(remainder, 3600)
    minutes, seconds = divrem(remainder, 60)

    if days > 0
        return @sprintf("%dd %02dh %02dm %02ds", days, hours, minutes, seconds)
    elseif hours > 0
        return @sprintf("%dh %02dm %02ds", hours, minutes, seconds)
    elseif minutes > 0
        return @sprintf("%dm %02ds", minutes, seconds)
    else
        return @sprintf("%ds", seconds)
    end
end



# ==============================================================================
# Helper: Dynamics Dispatcher
# ==============================================================================

function Step_Dynamics!(
    config,
    mesh,
    atmo_data,
    dyn_data,
    integrator,
    semi_implicit,
    vert_coord,
    physics_params,
)

    model_type = config.model_type
    if model_type == :Barotropic
        Barotropic_Dynamics!(mesh, atmo_data, dyn_data, integrator)

    elseif model_type == :Shallow_Water
        # Extract SW params safely
        kappa_m = get(physics_params, "kappa_m", 1.0 / (20.0 * 86400))
        kappa_t = get(physics_params, "kappa_t", 1.0 / (10.0 * 86400))
        h_eq = dyn_data.grid_geopots # We stored h_eq here during init
        h_0 = get(physics_params, "h_0", 3.0e4)

        Shallow_Water_Physics!(dyn_data, kappa_m, kappa_t, h_eq)
        Shallow_Water_Dynamics!(mesh, atmo_data, h_0, dyn_data, integrator)

    elseif model_type == :PrimitiveEquation
        # Primitive Equation Dynamics
        # Note: Atmosphere_Update! handles physics internals based on Atmo_Data flags
        Atmosphere_Update!(
            config,
            mesh,
            vert_coord,
            atmo_data,
            dyn_data,
            semi_implicit,
            physics_params,
        )

    else
        error("Unknown Model Type: $model_type")
    end
end

end
