using ...Experiment_Configuration
using ...Spectral_Spherical_Mesh_Module
using ...Vert_Coordinate_Module
using ...Atmo_Data_Module
using ...Dyn_Data_Module
using ...Semi_Implicit_Module
using ...Press_And_Geopot_Module

"""Persistent storage for the sequential, time-split gridpoint physics."""
struct Physics_Workspace
    grid_u::Array{Float64,3}
    grid_v::Array{Float64,3}
    grid_t::Array{Float64,3}
    grid_q::Array{Float64,3}
    grid_q_before::Array{Float64,3}
    grid_Δp_before::Array{Float64,3}
    grid_bm_t_tendency::Array{Float64,3}
    grid_bm_q_tendency::Array{Float64,3}
    grid_bm_precip::Array{Float64,3}
    grid_lscale_t_tendency::Array{Float64,3}
    grid_lscale_q_tendency::Array{Float64,3}
    grid_liquid_water_content::Array{Float64,3}
    grid_precip::Array{Float64,3}
    grid_shflx::Array{Float64,3}
    grid_lhflx::Array{Float64,3}
    grid_lrf_tendency::Array{Float64,3}
    pbl::PBL_Workspace
end

function Physics_Workspace(nλ::Int, nθ::Int, nd::Int)
    field() = zeros(Float64, nλ, nθ, nd)
    column() = zeros(Float64, nλ, nθ, 1)
    return Physics_Workspace(
        field(), field(), field(), field(), field(), field(),
        field(), field(), column(),
        field(), field(), field(), column(),
        column(), column(), field(),
        PBL_Workspace(nλ, nθ, nd),
    )
end

function _physics_workspace!(physics_params::Dict{String,Any}, atmo_data::Atmo_Data)
    workspace = get!(physics_params, "Physics_workspace") do
        Physics_Workspace(atmo_data.nλ, atmo_data.nθ, atmo_data.nd)
    end
    workspace isa Physics_Workspace ||
        error("physics_params[\"Physics_workspace\"] has an incompatible type")
    size(workspace.grid_t) == (atmo_data.nλ, atmo_data.nθ, atmo_data.nd) ||
        throw(DimensionMismatch("Physics_workspace does not match the model grid"))
    return workspace
end

function _reset_physics_diagnostics!(dyn_data::Dyn_Data)
    for field in (
        dyn_data.grid_shflx, dyn_data.grid_lhflx, dyn_data.grid_precip,
        dyn_data.grid_liquid_water_content, dyn_data.grid_t_eq,
        dyn_data.grid_lrf_tendency, dyn_data.grid_bm_t_tendency,
        dyn_data.grid_bm_q_tendency, dyn_data.grid_bm_precip,
    )
        fill!(field, 0.0)
    end
    return nothing
end

function _reset_substep_diagnostics!(workspace::Physics_Workspace)
    for field in (
        workspace.grid_bm_t_tendency, workspace.grid_bm_q_tendency,
        workspace.grid_bm_precip, workspace.grid_lscale_t_tendency,
        workspace.grid_lscale_q_tendency, workspace.grid_liquid_water_content,
        workspace.grid_precip, workspace.grid_shflx, workspace.grid_lhflx,
        workspace.grid_lrf_tendency,
    )
        fill!(field, 0.0)
    end
    return nothing
end

"""
    Spectral_Physics!(config, mesh, vert_coord, atmo_data, dyn_data,
                      semi_implicit, physics_params)

Apply physics directly to the provisional next state produced by the dynamics.
The leapfrog interval is covered by one or more `config.Δt` substeps, so no
physics increment is inserted into the leapfrog RHS or applied twice.

Within every substep the process order is convection, large-scale condensation,
surface exchange, PBL mixing, Rayleigh friction (with frictional heating),
Newtonian relaxation, and finally the moisture LRF.
"""
function Spectral_Physics!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    semi_implicit::Semi_Implicit_Solver,
    physics_params::Dict{String,Any},
)
    _reset_physics_diagnostics!(dyn_data)
    workspace = _physics_workspace!(physics_params, atmo_data)
    work_u, work_v, work_t, work_q =
        workspace.grid_u, workspace.grid_v, workspace.grid_t, workspace.grid_q
    copyto!(work_u, dyn_data.grid_u_n)
    copyto!(work_v, dyn_data.grid_v_n)
    copyto!(work_t, dyn_data.grid_t_n)
    copyto!(work_q, dyn_data.grid_q_n)

    effective_dt = Get_Δt(semi_implicit.integrator)
    physics_dt = config.Δt
    physics_dt > 0 || throw(ArgumentError("config.Δt must be positive"))
    effective_dt % physics_dt == 0 || throw(
        ArgumentError(
            "leapfrog interval $effective_dt s is not divisible by " *
            "physics timestep $physics_dt s",
        ),
    )
    nsubsteps = effective_dt ÷ physics_dt
    diagnostic_weight = 1.0 / nsubsteps

    grid_p_half = dyn_data.grid_p_half
    grid_p_full = dyn_data.grid_p_full
    grid_ps = dyn_data.grid_ps_n
    do_betts_miller = get(physics_params, "do_Betts_Miller", false)
    do_lscale_cond = get(physics_params, "do_Lscale_Cond", false)
    do_sensible = get(physics_params, "do_Sensible_Heating", false)
    do_evaporation =
        get(physics_params, "do_Surface_Evaporation", false) &&
        config.moisture_processes
    do_pbl_mixing = get(physics_params, "do_Implicit_PBL_Scheme", false)
    do_hs = get(physics_params, "do_HS_Forcing", false)
    do_lrf = get(physics_params, "do_LRF", false)

    if do_betts_miller || do_lscale_cond
        config.moisture_processes ||
            error("moist convection and condensation require moisture_processes = true")
    end
    if do_betts_miller
        haskey(physics_params, "BM_state") ||
            error("BM_state was not initialized before time integration")
        bm_state = physics_params["BM_state"]::Betts_Miller_State
        physics_dt <= bm_state.tau || throw(
            ArgumentError(
                "Betts-Miller requires physics_dt <= bm_tau; got $physics_dt s and $(bm_state.tau) s",
            ),
        )
    else
        bm_state = nothing
    end
    if do_lrf
        config.moisture_processes || error("LRF requires moisture_processes = true")
        haskey(physics_params, "LRF_state") ||
            error("LRF_state was not initialized before time integration")
    end

    heating_fraction = get(
        physics_params,
        "condensation_heating_fraction",
        get(physics_params, "L", 1.0),
    )
    lower_boundary_temperature = get(
        physics_params,
        "lower_boundary_temperature",
        Default_Lower_Boundary_Temperature,
    )

    for _ in 1:nsubsteps
        _reset_substep_diagnostics!(workspace)
        copyto!(workspace.grid_q_before, work_q)
        @inbounds for k in 1:atmo_data.nd
            @views @. workspace.grid_Δp_before[:, :, k] =
                grid_p_half[:, :, k+1] - grid_p_half[:, :, k]
        end

        if do_betts_miller || do_lscale_cond
            Moist_Physics!(
                do_betts_miller, do_lscale_cond, bm_state, heating_fraction,
                atmo_data, work_t, work_q, grid_p_full, grid_p_half, physics_dt,
                workspace.grid_bm_t_tendency, workspace.grid_bm_q_tendency,
                workspace.grid_bm_precip, workspace.grid_lscale_t_tendency,
                workspace.grid_lscale_q_tendency, workspace.grid_liquid_water_content,
                workspace.grid_precip,
            )
            @. dyn_data.grid_bm_t_tendency +=
                diagnostic_weight * workspace.grid_bm_t_tendency
            @. dyn_data.grid_bm_q_tendency +=
                diagnostic_weight * workspace.grid_bm_q_tendency
            @. dyn_data.grid_bm_precip += diagnostic_weight * workspace.grid_bm_precip
            @. dyn_data.grid_liquid_water_content +=
                diagnostic_weight * workspace.grid_liquid_water_content
            @. dyn_data.grid_precip += diagnostic_weight * workspace.grid_precip
        end

        if do_sensible || do_evaporation || do_pbl_mixing
            surface_wind, surface_height, density = Calculate_V_c_za_rho!(
                workspace.pbl, atmo_data, grid_p_half, grid_p_full, grid_ps,
                work_u, work_v, work_t, work_q,
            )
            if do_sensible
                Sensible_Heating!(
                    mesh, atmo_data, work_t, workspace.grid_shflx, surface_wind,
                    surface_height, density, physics_dt,
                    physics_params["C_H"]::Float64, lower_boundary_temperature,
                )
                @. dyn_data.grid_shflx += diagnostic_weight * workspace.grid_shflx
            end
            if do_evaporation
                Surface_Evaporation!(
                    mesh, atmo_data, grid_ps, work_q, workspace.grid_lhflx,
                    surface_wind, surface_height, density, physics_dt,
                    physics_params["C_E"]::Float64, lower_boundary_temperature,
                )
                @. dyn_data.grid_lhflx += diagnostic_weight * workspace.grid_lhflx
            end
            if do_pbl_mixing
                Implicit_PBL_Mixing!(
                    atmo_data, grid_p_full, grid_p_half, work_t, work_q,
                    dyn_data.K_E, surface_wind, surface_height, physics_params,
                    physics_dt, physics_params["C_D"]::Float64,
                )
            end
        else
            fill!(dyn_data.K_E, 0.0)
        end

        if do_hs
            Rayleigh_Friction!(
                atmo_data, physics_dt, config.day_to_sec, grid_p_half, grid_p_full,
                work_u, work_v, work_t, physics_params,
            )
            Newtonian_Relaxation!(
                atmo_data, physics_dt, config.day_to_sec, mesh.sinθ, grid_p_half,
                grid_p_full, work_t, dyn_data.grid_t_eq, physics_params,
            )
        end

        if do_lrf
            LRF!(
                physics_params["LRF_state"]::LRF_State,
                work_q, workspace.grid_lrf_tendency, config.day_to_sec,
            )
            @. work_t += physics_dt * workspace.grid_lrf_tendency
            @. dyn_data.grid_lrf_tendency +=
                diagnostic_weight * workspace.grid_lrf_tendency
        end

        all(isfinite, work_u) && all(isfinite, work_v) && all(isfinite, work_t) &&
            all(isfinite, work_q) ||
            error("physics produced a non-finite prognostic state")
        minimum(work_t) > 0.0 || error("physics produced a non-positive temperature")
        minimum(work_q) >= -1.0e-14 || error("physics produced negative specific humidity")
        work_q .= max.(work_q, 0.0)
        maximum(work_q) < 1.0 || error("physics produced specific humidity >= 1")

        if config.moisture_processes
            Dry_Air_Adjustment!(
                vert_coord,
                grid_ps,
                workspace.grid_q_before,
                work_q,
                workspace.grid_Δp_before,
            )
            Pressure_Variables!(
                vert_coord,
                grid_ps,
                dyn_data.grid_p_half,
                dyn_data.grid_Δp,
                dyn_data.grid_lnp_half,
                dyn_data.grid_p_full,
                dyn_data.grid_lnp_full,
            )
        end
    end

    copyto!(dyn_data.grid_u_n, work_u)
    copyto!(dyn_data.grid_v_n, work_v)
    copyto!(dyn_data.grid_t_n, work_t)
    copyto!(dyn_data.grid_q_n, work_q)
    return nothing
end
