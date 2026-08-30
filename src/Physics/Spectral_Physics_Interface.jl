using ...Experiment_Configuration
using ...Spectral_Spherical_Mesh_Module
using ...Vert_Coordinate_Module
using ...Atmo_Data_Module
using ...Dyn_Data_Module
using ...Semi_Implicit_Module

"""Persistent storage for the sequential gridpoint physics calculation."""
struct Physics_Workspace
    grid_u::Array{Float64,3}
    grid_v::Array{Float64,3}
    grid_t::Array{Float64,3}
    grid_q::Array{Float64,3}
    grid_lscale_t_tendency::Array{Float64,3}
    grid_lscale_q_tendency::Array{Float64,3}
    pbl::PBL_Workspace
end

function Physics_Workspace(nλ::Int, nθ::Int, nd::Int)
    field() = zeros(Float64, nλ, nθ, nd)
    return Physics_Workspace(
        field(),
        field(),
        field(),
        field(),
        field(),
        field(),
        PBL_Workspace(nλ, nθ, nd),
    )
end

function _physics_workspace!(
    physics_params::Dict{String,Any},
    atmo_data::Atmo_Data,
)
    workspace = get!(physics_params, "Physics_workspace") do
        Physics_Workspace(atmo_data.nλ, atmo_data.nθ, atmo_data.nd)
    end
    workspace isa Physics_Workspace ||
        error("physics_params[\"Physics_workspace\"] has an incompatible type")
    size(workspace.grid_t) == (atmo_data.nλ, atmo_data.nθ, atmo_data.nd) ||
        throw(DimensionMismatch("Physics_workspace does not match the model grid"))
    return workspace
end

"""
    Spectral_Physics!(config, mesh, vert_coord, atmo_data, dyn_data,
                      semi_implicit, physics_params)

Apply all physical parameterizations sequentially to a private working state,
then expose one combined process-split tendency to the leapfrog dynamical core.
The prognostic current grid and spectral states are never modified here.
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
    grid_δu = dyn_data.grid_δu
    grid_δv = dyn_data.grid_δv
    grid_δps = dyn_data.grid_δps
    grid_δt = dyn_data.grid_δt
    grid_δq = dyn_data.grid_δq

    # All tendency and per-step diagnostic arrays are overwritten once here.
    fill!(grid_δu, 0.0)
    fill!(grid_δv, 0.0)
    fill!(grid_δps, 0.0)
    fill!(grid_δt, 0.0)
    fill!(grid_δq, 0.0)
    fill!(dyn_data.spe_δq, 0.0)
    fill!(dyn_data.grid_shflx, 0.0)
    fill!(dyn_data.grid_lhflx, 0.0)
    fill!(dyn_data.grid_precip, 0.0)
    fill!(dyn_data.grid_liquid_water_content, 0.0)
    fill!(dyn_data.grid_t_eq, 0.0)
    fill!(dyn_data.grid_lrf_tendency, 0.0)
    fill!(dyn_data.grid_bm_t_tendency, 0.0)
    fill!(dyn_data.grid_bm_q_tendency, 0.0)
    fill!(dyn_data.grid_bm_precip, 0.0)

    workspace = _physics_workspace!(physics_params, atmo_data)
    work_u = workspace.grid_u
    work_v = workspace.grid_v
    work_t = workspace.grid_t
    work_q = workspace.grid_q
    copyto!(work_u, dyn_data.grid_u_c)
    copyto!(work_v, dyn_data.grid_v_c)
    copyto!(work_t, dyn_data.grid_t_c)
    copyto!(work_q, dyn_data.grid_q_c)

    effective_dt = Get_Δt(semi_implicit.integrator)
    grid_p_half = dyn_data.grid_p_half
    grid_p_full = dyn_data.grid_p_full
    grid_ps = dyn_data.grid_ps_c

    # Convection is followed by a grid-scale saturation cleanup. Both schemes
    # update only the private working state and retain separate diagnostics.
    do_betts_miller = get(physics_params, "do_Betts_Miller", false)
    do_lscale_cond = get(physics_params, "do_Lscale_Cond", false)
    if do_betts_miller || do_lscale_cond
        config.moisture_processes ||
            error("moist convection and condensation require moisture_processes = true")

        if do_betts_miller
            haskey(physics_params, "BM_state") ||
                error("BM_state was not initialized before time integration")
            bm_state = physics_params["BM_state"]::Betts_Miller_State
            effective_dt <= bm_state.tau || throw(
                ArgumentError(
                    "Betts-Miller requires effective_dt <= bm_tau; " *
                    "got $effective_dt s and $(bm_state.tau) s",
                ),
            )
        else
            bm_state = nothing
        end

        heating_fraction = get(
            physics_params,
            "condensation_heating_fraction",
            get(physics_params, "L", 1.0),
        )
        Moist_Physics!(
            do_betts_miller,
            do_lscale_cond,
            bm_state,
            heating_fraction,
            atmo_data,
            work_t,
            work_q,
            grid_p_full,
            grid_p_half,
            effective_dt,
            dyn_data.grid_bm_t_tendency,
            dyn_data.grid_bm_q_tendency,
            dyn_data.grid_bm_precip,
            workspace.grid_lscale_t_tendency,
            workspace.grid_lscale_q_tendency,
            dyn_data.grid_liquid_water_content,
            dyn_data.grid_precip,
        )
    else
        fill!(workspace.grid_lscale_t_tendency, 0.0)
        fill!(workspace.grid_lscale_q_tendency, 0.0)
    end

    # Surface exchange and vertical mixing see the post-condensation state.
    do_sensible = get(physics_params, "do_Sensible_Heating", false)
    do_evaporation =
        get(physics_params, "do_Surface_Evaporation", false) && config.moisture_processes
    do_pbl_mixing = get(physics_params, "do_Implicit_PBL_Scheme", false)
    if do_sensible || do_evaporation || do_pbl_mixing
        surface_wind, surface_height, density = Calculate_V_c_za_rho!(
            workspace.pbl,
            atmo_data,
            grid_p_half,
            grid_p_full,
            grid_ps,
            work_u,
            work_v,
            work_t,
            work_q,
        )
        lower_boundary_temperature = get(
            physics_params,
            "lower_boundary_temperature",
            Default_Lower_Boundary_Temperature,
        )

        if do_sensible
            Sensible_Heating!(
                mesh,
                atmo_data,
                work_t,
                dyn_data.grid_shflx,
                surface_wind,
                surface_height,
                density,
                effective_dt,
                physics_params["C_H"]::Float64,
                lower_boundary_temperature,
            )
        end
        if do_evaporation
            Surface_Evaporation!(
                mesh,
                atmo_data,
                grid_ps,
                work_q,
                dyn_data.grid_lhflx,
                surface_wind,
                surface_height,
                density,
                effective_dt,
                physics_params["C_E"]::Float64,
                lower_boundary_temperature,
            )
        end
        if do_pbl_mixing
            Implicit_PBL_Mixing!(
                atmo_data,
                grid_p_full,
                grid_p_half,
                work_t,
                work_q,
                dyn_data.K_E,
                surface_wind,
                surface_height,
                physics_params,
                effective_dt,
                physics_params["C_D"]::Float64,
            )
        end
    else
        fill!(dyn_data.K_E, 0.0)
    end

    # Rayleigh damping, including exact finite-step frictional heating, precedes
    # the Newtonian relaxation so radiation sees the post-friction temperature.
    if get(physics_params, "do_HS_Forcing", false)
        Rayleigh_Friction!(
            atmo_data,
            effective_dt,
            config.day_to_sec,
            grid_p_half,
            grid_p_full,
            work_u,
            work_v,
            work_t,
            physics_params,
        )
        Newtonian_Relaxation!(
            atmo_data,
            effective_dt,
            config.day_to_sec,
            mesh.sinθ,
            grid_p_half,
            grid_p_full,
            work_t,
            dyn_data.grid_t_eq,
            physics_params,
        )
    end

    # The moisture LRF is the final radiative contribution and therefore sees
    # the humidity left by convection, condensation, surface exchange and PBL.
    if get(physics_params, "do_LRF", false)
        config.moisture_processes || error("LRF requires moisture_processes = true")
        haskey(physics_params, "LRF_state") ||
            error("LRF_state was not initialized before time integration")
        LRF!(
            physics_params["LRF_state"]::LRF_State,
            work_q,
            dyn_data.grid_lrf_tendency,
            config.day_to_sec,
        )
        @. work_t += effective_dt * dyn_data.grid_lrf_tendency
    end

    # This telescoping state difference is the sole authoritative physics
    # tendency. It includes each sequential process exactly once.
    inv_effective_dt = 1.0 / Float64(effective_dt)
    @. grid_δu = (work_u - dyn_data.grid_u_c) * inv_effective_dt
    @. grid_δv = (work_v - dyn_data.grid_v_c) * inv_effective_dt
    @. grid_δt = (work_t - dyn_data.grid_t_c) * inv_effective_dt
    if config.moisture_processes
        @. grid_δq = (work_q - dyn_data.grid_q_c) * inv_effective_dt
    end

    return nothing
end
