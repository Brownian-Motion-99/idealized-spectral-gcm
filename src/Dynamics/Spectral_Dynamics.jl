module Spectral_Dynamics_Module

using Base.Threads
using LinearAlgebra
using Printf
using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Vert_Coordinate_Module
using ..Dyn_Data_Module
using ..Semi_Implicit_Module
using ..Time_Integrator_Module
using ..Grid_Tracer_Transport_Module
using ..Press_And_Geopot_Module
using ..Atmos_Param_Module
using ..Experiment_Configuration

using Statistics
using Interpolations
export Compute_Corrections_Init,
    Compute_Corrections!,
    Reset_Dynamics_Tendencies!,
    Four_In_One!,
    Spectral_Dynamics!,
    Get_Topography!,
    Spectral_Initialize_Fields!,
    Spectral_Dynamics_Physics!,
    Atmosphere_Update!



"""
    Compute_Corrections_Init(
        mesh, vert_coord, atmo_data,
        grid_ps_p,
        grid_energy_temp,
        grid_u, grid_v, grid_t, grid_q
    )

Calculates the global integrals of mass, total energy, and moisture to serve as 
conservation targets for a single, internally consistent model state. Physics is
applied after the dynamics update, so physics tendencies must not be folded into
these dynamics-only targets.

### Parameters
    - mesh: Spectral spherical mesh properties.
    - vert_coord: Vertical coordinate definitions.
    - atmo_data: Atmospheric physical constants and flags (do_mass_correction, etc.).

    - grid_ps_p: Surface pressure of the target state [nλ, nθ, 1].
    
    - grid_energy_temp: Temporary work array for energy calculation [nλ, nθ, nd].
    
    - grid_u, grid_v: Horizontal wind components of the target state.
    - grid_t: Temperature of the target state.
    - grid_q: Specific humidity of the target state.

### Returns
    - mean_ps_p: Global mean surface pressure (Mass target).
    - mean_energy_p: Global integral of total energy (Enthalpy + Kinetic).
    - mean_moisture_p: Global integral of specific humidity.

"""
function Compute_Corrections_Init(
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    grid_ps_p::Array{Float64,3},
    grid_energy_temp::Array{Float64,3},
    grid_u::Array{Float64,3},
    grid_v::Array{Float64,3},
    grid_t::Array{Float64,3},
    grid_q::Array{Float64,3},
)

    mean_ps_p = 0.0
    mean_energy_p = 0.0
    mean_moisture_p = 0.0

    do_mass_correction, do_energy_correction, do_water_correction =
        atmo_data.do_mass_correction,
        atmo_data.do_energy_correction,
        atmo_data.do_water_correction
    cp_air = atmo_data.cp_air

    if (do_mass_correction)
        mean_ps_p = Area_Weighted_Global_Mean(mesh, grid_ps_p)
    end

    if (do_energy_correction)
        @. grid_energy_temp = 0.5 * (grid_u^2 + grid_v^2) + cp_air * grid_t
        mean_energy_p = Mass_Weighted_Global_Integral(
            vert_coord,
            mesh,
            atmo_data,
            grid_energy_temp,
            grid_ps_p,
        )
    end

    if (do_water_correction)
        mean_moisture_p = Mass_Weighted_Global_Integral(
            vert_coord,
            mesh,
            atmo_data,
            grid_q,
            grid_ps_p,
        )
    end

    return mean_ps_p, mean_energy_p, mean_moisture_p
end

"""Clear every explicit dynamics tendency before constructing a new RHS."""
function Reset_Dynamics_Tendencies!(dyn_data::Dyn_Data)
    for field in (
        dyn_data.spe_δvor,
        dyn_data.spe_δdiv,
        dyn_data.spe_δlnps,
        dyn_data.spe_δt,
        dyn_data.grid_δvor,
        dyn_data.grid_δdiv,
        dyn_data.grid_δu,
        dyn_data.grid_δv,
        dyn_data.grid_δps,
        dyn_data.grid_δlnps,
        dyn_data.grid_δt,
        dyn_data.grid_δq,
    )
        fill!(field, zero(eltype(field)))
    end
    return nothing
end

"""Apply the Robert filter to grid humidity without changing its current-state water."""
function _filter_grid_humidity!(
    integrator::Filtered_Leapfrog,
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    grid_q_p::Array{Float64,3},
    grid_q_c::Array{Float64,3},
    grid_q_n::Array{Float64,3},
    grid_ps_c::Array{Float64,3},
)
    target_water = atmo_data.do_water_correction ?
        Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_q_c, grid_ps_c) : 0.0
    α = integrator.robert_coef
    if integrator.init_step
        @. grid_q_c += α * (-grid_q_c + grid_q_n)
    else
        @. grid_q_c += α * (grid_q_p - 2grid_q_c + grid_q_n)
    end
    if atmo_data.do_water_correction
        Restore_Grid_Tracer_Integral!(
            grid_q_c, grid_ps_c, target_water, vert_coord, mesh, atmo_data,
        )
    end
    return nothing
end

"""Make a gridpoint-adjusted next state authoritative in spectral space."""
function _synchronize_physics_next!(
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
)
    mean_ps, mean_energy, _ = Compute_Corrections_Init(
        mesh,
        vert_coord,
        atmo_data,
        dyn_data.grid_ps_n,
        dyn_data.grid_energy_full,
        dyn_data.grid_u_n,
        dyn_data.grid_v_n,
        dyn_data.grid_t_n,
        dyn_data.grid_q_n,
    )
    mean_water = Mass_Weighted_Global_Integral(
        vert_coord, mesh, atmo_data, dyn_data.grid_q_n, dyn_data.grid_ps_n,
    )

    all(isfinite, dyn_data.grid_ps_n) && minimum(dyn_data.grid_ps_n) > 0.0 ||
        error("physics produced invalid surface pressure")
    @. dyn_data.grid_lnps = log(dyn_data.grid_ps_n)
    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_lnps, dyn_data.spe_lnps_n)
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_lnps_n, dyn_data.grid_lnps)
    @. dyn_data.grid_ps_n = exp(dyn_data.grid_lnps)

    Vor_Div_From_Grid_UV!(
        mesh,
        dyn_data.grid_u_n,
        dyn_data.grid_v_n,
        dyn_data.spe_vor_n,
        dyn_data.spe_div_n,
    )
    UV_Grid_From_Vor_Div!(
        mesh,
        dyn_data.spe_vor_n,
        dyn_data.spe_div_n,
        dyn_data.grid_u_n,
        dyn_data.grid_v_n,
    )
    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_t_n, dyn_data.spe_t_n)
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_t_n, dyn_data.grid_t_n)
    Compute_Corrections!(
        mesh,
        vert_coord,
        atmo_data,
        mean_ps,
        dyn_data.grid_ps_n,
        dyn_data.spe_lnps_n,
        mean_energy,
        dyn_data.grid_energy_full,
        dyn_data.grid_u_n,
        dyn_data.grid_v_n,
        dyn_data.grid_t_n,
        dyn_data.spe_t_n,
        mean_water,
        dyn_data.grid_q_n,
        ;
        correct_water = false,
    )

    # Pressure projection and the mass/energy fixers are now complete. Humidity
    # is grid-only, and this must be its last synchronization operation.
    Restore_Grid_Tracer_Integral!(
        dyn_data.grid_q_n,
        dyn_data.grid_ps_n,
        mean_water,
        vert_coord,
        mesh,
        atmo_data,
    )

    all(isfinite, dyn_data.grid_u_n) && all(isfinite, dyn_data.grid_v_n) &&
        all(isfinite, dyn_data.grid_t_n) && all(isfinite, dyn_data.grid_q_n) ||
        error("spectral projection produced a non-finite post-physics state")
    minimum(dyn_data.grid_t_n) > 0.0 ||
        error("spectral projection produced a non-positive temperature")
    minimum(dyn_data.grid_q_n) >= -1.0e-14 ||
        error("post-physics synchronization produced negative specific humidity")
    return nothing
end



"""
    Compute_Corrections!(
        mesh, vert_coord, atmo_data,
        mean_ps_p, grid_ps_n, spe_lnps_n, 
        mean_energy_p, grid_energy_temp, 
        grid_u_n, grid_v_n, grid_t_n, spe_t_n,
        mean_moisture_p, grid_q_n
    )

Adjusts the updated atmospheric state fields (grid_*_n) to strictly enforce global 
conservation laws for mass, energy, and water vapor. The correction is applied 
multiplicatively for mass and moisture, and additively for temperature (energy).

### Parameters
    - mesh: Spectral spherical mesh properties.
    - vert_coord: Vertical coordinate definitions.
    - atmo_data: Atmospheric physical constants.

    - mean_ps_p: Target global mean surface pressure (computed in `Compute_Corrections_Init`).
    - grid_ps_n: Updated surface pressure field at (t+Δt). Modified in-place.
    - spe_lnps_n: Spectral coefficient of log-surface pressure (0,0 mode modified).
    
    - mean_energy_p: Target global total energy integral.
    - grid_energy_temp: Temporary work array.
    
    - grid_u_n, grid_v_n: Updated horizontal wind components at (t+Δt).
    - grid_t_n: Updated temperature field at (t+Δt). Modified in-place.
    - spe_t_n: Spectral coefficients of temperature. Modified in-place.
    
    - mean_moisture_p: Target global moisture integral.
    - grid_q_n: Updated specific humidity field at (t+Δt). Modified in-place.

### Returns
    - nothing

### Modified
    - grid_ps_n
    - spe_lnps_n
    - grid_t_n
    - spe_t_n
    - grid_q_n

"""
function Compute_Corrections!(
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    mean_ps_p::Float64,
    grid_ps_n::Array{Float64,3},
    spe_lnps_n::Array{ComplexF64,3},
    mean_energy_p::Float64,
    grid_energy_temp::Array{Float64,3},
    grid_u_n::Array{Float64,3},
    grid_v_n::Array{Float64,3},
    grid_t_n::Array{Float64,3},
    spe_t_n::Array{ComplexF64,3},
    mean_moisture_p::Float64,
    grid_q_n::Array{Float64,3},
    ;
    correct_water::Bool = true,
)

    do_mass_correction, do_energy_correction, do_water_correction =
        atmo_data.do_mass_correction,
        atmo_data.do_energy_correction,
        atmo_data.do_water_correction
    cp_air, grav = atmo_data.cp_air, atmo_data.grav

    if (do_mass_correction)
        mean_ps_n = Area_Weighted_Global_Mean(mesh, grid_ps_n)
        mass_correction_factor = mean_ps_p / mean_ps_n
        grid_ps_n .*= mass_correction_factor
        spe_lnps_n[1, 1, 1] += log(mass_correction_factor)
    end

    if (do_energy_correction)
        grid_energy_temp .= 0.5 * (grid_u_n .^ 2 + grid_v_n .^ 2) + cp_air * grid_t_n
        mean_energy_n = Mass_Weighted_Global_Integral(
            vert_coord,
            mesh,
            atmo_data,
            grid_energy_temp,
            grid_ps_n,
        )

        mean_ps_for_energy = Area_Weighted_Global_Mean(mesh, grid_ps_n)
        isfinite(mean_ps_for_energy) && mean_ps_for_energy > 0.0 ||
            error("Cannot apply energy correction: surface pressure is non-positive")
        temperature_correction =
            grav * (mean_energy_p - mean_energy_n) / (cp_air * mean_ps_for_energy)
        grid_t_n .+= temperature_correction
        spe_t_n[1, 1, :] .+= temperature_correction
    end

    if do_water_correction && correct_water
        Restore_Grid_Tracer_Integral!(
            grid_q_n, grid_ps_n, mean_moisture_p, vert_coord, mesh, atmo_data,
        )
    end

end



"""
    Four_In_One!(
        vert_coord, atmo_data, 
        grid_div, grid_u, grid_v, 
        grid_ps, grid_Δp, grid_lnp_half, grid_lnp_full, grid_p_full,
        grid_dλ_ps, grid_dθ_ps, 
        grid_virtual_t,
        grid_M_half, grid_w_full, 
        grid_δu, grid_δv, grid_δps, grid_δt
    )

Simultaneously computes the vertical mass flux (M), pressure vertical velocity (ω), 
and updates the tendencies for temperature, momentum, and surface pressure based 
on the adiabatic primitive equations.

### Parameters
    - vert_coord: Vertical discretization parameters.
    - atmo_data: Physical constants.

    - grid_div: Divergence field [nλ, nθ, nd].
    - grid_u, grid_v: Zonal and meridional wind components.
    
    - grid_ps: Surface pressure.
    - grid_Δp: Pressure thickness of layers.
    - grid_lnp_half, grid_lnp_full: Log-pressure at interfaces and layer centers.
    - grid_p_full: Pressure at layer centers.
    
    - grid_dλ_ps, grid_dθ_ps: Zonal and meridional gradients of surface pressure.
    
    - grid_virtual_t: Equation-of-state temperature field.
    
    - grid_M_half: Output vertical mass flux at interfaces [nλ, nθ, nd+1].
    - grid_w_full: Output vertical velocity (ω = dp/dt) at layer centers [nλ, nθ, nd].
    
    - grid_δu, grid_δv: Momentum tendencies. Modified in-place (pressure gradient force subtracted).
    - grid_δps: Surface pressure tendency. Modified in-place (continuity equation).
    - grid_δt: Temperature tendency. Modified in-place (energy conversion added).
### Returns
    - nothing

### Modified
    - grid_M_half
    - grid_w_full
    - grid_δu, grid_δv
    - grid_δps
    - grid_δt

### Notes
    - **Continuity Equation (Surface Pressure):**
    ∂pₛ/∂t = -∑ₖ ∇⋅(𝐯ₖ Δpₖ)
    The function accumulates the divergence of mass flux (`dmean`) down the column 
    to find the total surface pressure tendency (`dmean_tot`).

    - **Vertical Mass Flux (M):**
    Represents the vertical mass flux relative to the hybrid coordinate surfaces:
    Mₖ₊½ = -∑_{r=1}ᵏ ∇⋅(𝐯ᵣ Δpᵣ) - Bₖ₊½ (∂pₛ/∂t)

    - **Pressure Gradient Force (Momentum):**
    Updates the momentum tendencies by subtracting the gradient of geopotential/pressure:
    δu -= RT ⋅ ∇lnp
    The gradient ∇lnp is computed locally using the chain rule on surface pressure 
    gradients and vertical coordinate coefficients.

    - **Energy Conversion (Temperature):**
    Updates the temperature tendency with the adiabatic expansion/compression term:
    δT += κ T ω / p
    
    - **Vertical Velocity (ω):**
    Computed as the total derivative of pressure:
    ω = dp/dt = ∂p/∂t + 𝐯⋅∇p + M ∂p/∂σ

"""
function Four_In_One!(
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    grid_div::Array{Float64,3},
    grid_u::Array{Float64,3},
    grid_v::Array{Float64,3},
    grid_ps::Array{Float64,3},
    grid_Δp::Array{Float64,3},
    grid_lnp_half::Array{Float64,3},
    grid_lnp_full::Array{Float64,3},
    grid_p_full::Array{Float64,3},
    grid_dλ_ps::Array{Float64,3},
    grid_dθ_ps::Array{Float64,3},
    grid_virtual_t::Array{Float64,3},
    grid_M_half::Array{Float64,3},
    grid_w_full::Array{Float64,3},
    grid_δu::Array{Float64,3},
    grid_δv::Array{Float64,3},
    grid_δps::Array{Float64,3},
    grid_δt::Array{Float64,3},
)

    # Unpack parameters
    rdgas, cp_air = atmo_data.rdgas, atmo_data.cp_air
    nd, bk = vert_coord.nd, vert_coord.bk
    Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
    vert_difference_option = vert_coord.vert_difference_option
    kappa = rdgas / cp_air
    nλ, nθ = atmo_data.nλ, atmo_data.nθ

    @threads for j = 1:nθ
        for i = 1:nλ

            if (vert_difference_option == "simmons_and_burridge")

                # dmean_tot = ∇ ∑_{k=1}^{nd} v_k * Δp_k = ∑_{k=1}^{nd} D_k
                dmean_tot = 0.0

                for k = 1:nd

                    Δp = grid_Δp[i, j, k]
                    Δlnp_p = grid_lnp_half[i, j, k+1] - grid_lnp_full[i, j, k]
                    Δlnp_m = grid_lnp_full[i, j, k] - grid_lnp_half[i, j, k]
                    Δlnp = grid_lnp_half[i, j, k+1] - grid_lnp_half[i, j, k]

                    # Pressure gradient force
                    # ∇p_k/p = [(lnp_k - lnp_{k-1/2}) * ∇p_{k-1/2} + (lnp_{k+1/2} - lnp_k) * ∇p_{k+1/2}] / Δpk
                    #        = [(lnp_k - lnp_{k-1/2}) * b_{k-1/2} + (lnp_{k+1/2} - lnp_k) * b_{k+1/2}] / Δpk * ∇ps
                    #        = x1 * ∇ps
                    x1 = (bk[k] * Δlnp_m + bk[k+1] * Δlnp_p) / Δp    # Pressure gradient mapper
                    dlnp_dλ = x1 * grid_dλ_ps[i, j, 1]
                    dlnp_dθ = x1 * grid_dθ_ps[i, j, 1]

                    # (grid_δu, grid_δv) -= R*T_v/p * ∇p
                    virtual_t_k = grid_virtual_t[i, j, k]
                    grid_δu[i, j, k] -= rdgas * virtual_t_k * dlnp_dλ
                    grid_δv[i, j, k] -= rdgas * virtual_t_k * dlnp_dθ

                    # dmean = ∇⋅(v_k * Δp_k) = div_k * Δp_k + v_k * (Δbk_k * ∇p_s)
                    dmean =
                        grid_div[i, j, k] * Δp +
                        Δbk[k] * (
                            grid_u[i, j, k] * grid_dλ_ps[i, j, 1] +
                            grid_v[i, j, k] * grid_dθ_ps[i, j, 1]
                        )

                    # energy conservation for temperature
                    # w/p = dlnp/dt = ∂lnp/∂t + ∂lnp/∂σ dσ + v⋅∇lnp
                    # ∂ξ_k/∂σ dσ = [M_{k+1/2}(ξ_k+1/2 - ξ_k) + M_{k-1/2}(ξ_k - ξ_k-1/2)]/Δp_k, where ξ is a dummy variable
                    # weight the same way
                    # vertical advection operator (M is the downward speed)
                    # ∂lnp_k/∂σ dσ = [M_{k+1/2} (lnp_k+1/2 - lnp_k) + M_{k-1/2} (lnp_k - lnp_k-1/2)] / Δp_k
                    # ∂lnp/∂t = 1/p ∂p/∂t = [∂p/∂t_{k+1/2} (lnp_k+1/2 - lnp_k) + ∂p/∂t_{k-1/2} (lnp_k - lnp_k-1/2)] / Δp_k
                    # As we know
                    # ∂p/∂t_{k+1/2} = -∑_{r=1}^k Dr - M_{k+1/2}
                    # Therefore
                    # ∂lnp/∂t + dσ ∂lnp/∂σ =  [(-∑_{r=1}^k Dr)(lnp_k+1/2 - lnp_k) + (-∑_{r=1}^{k-1} Dr)(lnp_k - lnp_k-1/2)]/Δp_k
                    #                      = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k
                    # w/p = x5 = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k + v⋅∇lnp
                    x5 =
                        -(dmean_tot * Δlnp + dmean * Δlnp_p) / Δp +
                        grid_u[i, j, k] * dlnp_dλ +
                        grid_v[i, j, k] * dlnp_dθ

                    # Constant-cp moist pressure work: grid_δt += κ*T_v*w/p
                    grid_δt[i, j, k] += kappa * virtual_t_k * x5

                    # grid_w_full[:,:,k] = dp/dt vertical velocity 
                    grid_w_full[i, j, k] = x5 * grid_p_full[i, j, k]

                    # update dmean_tot to ∑_{r=1}^k ∇⋅(v_r * Δp_r)
                    dmean_tot += dmean

                    # grid_M_half[:,:,k+1] = downward mass flux/per unit area across the K+1/2
                    # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
                    grid_M_half[i, j, k+1] = -dmean_tot

                end

            else
                error(
                    "vert_difference_option ",
                    vert_difference_option,
                    " is not a valid value for option",
                )
            end

            # ∂ps/∂t = -∑_{r=1}^nd ∇(vrΔp_r) = -dmean_tot
            grid_δps[i, j, 1] -= dmean_tot

            for k = 1:nd-1
                # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
                grid_M_half[i, j, k+1] += dmean_tot * bk[k+1]
            end

            grid_M_half[i, j, 1] = 0.0
            grid_M_half[i, j, nd+1] = 0.0

        end
    end

end



"""
    Spectral_Dynamics!(
        config,
        mesh, vert_coord, 
        atmo_data, dyn_data, 
        semi_implicit
    )

Performs one full time step of the primitive equations using a spectral transform method 
combined with a semi-implicit semi-Lagrangian (or Eulerian) time integration scheme.

### Parameters
    - config: Model configuration struct (physics switches, options).
    - mesh: Spectral spherical mesh properties (transforms, wavenumbers).
    - vert_coord: Vertical grid coordinates (hybrid sigma-pressure).
    - atmo_data: Atmospheric constants and parameters.
    - dyn_data: The main data container holding state variables for previous (p), current (c), and next (n) time steps.
    - semi_implicit: Solver structure for the implicit gravity wave correction.

### Returns
    - nothing

### Modified
    - dyn_data

### Notes
    - **Algorithm Sequence:**
    1. **Conservation Check:** Computes global mass/energy integrals of the previous state.
    2. **Pressure-related:** Updates pressure variables and computes gradients ∇pₛ.
    3. **Adiabatic Tendencies:** Calls `Four_In_One!` to compute vertical velocity,
       the moist pressure-gradient force, and constant-cp pressure work.
    4. **Vertical Advection:** Explicitly computes vertical advection for u, v, T, and tracers.
    5. **Horizontal Advection:** Computes non-linear advection terms.
    6. **Vector Invariant Formulation:** Computes vorticity and divergence tendencies (ηv, ηu) and kinetic energy (E = Φ + ½(u² + v²)).
    7. **Implicit Correction:** Inverts the Helmholtz equation in spectral space to stabilize gravity waves (`Implicit_Correction!`).
    8. **Spectral Damping:** Applies ∇²ⁿ hyper-diffusion to dampen small-scale noise.
    9. **Time Integration:** Advances the state using `Filtered_Leapfrog!` (Robert-Asselin-Williams filter).
    10. **Conservation Fixer:** Corrects the final predicted state to match global integrals.
    11. **Time Advance:** Rotates time indices (p ← c, c ← n).

"""
function Spectral_Dynamics!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    semi_implicit::Semi_Implicit_Solver,
    ;
    advance_time::Bool = true,
)

    Reset_Dynamics_Tendencies!(dyn_data)

    # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps =
        dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor =
        dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv =
        dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt =
        dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    # grid quantities
    grid_u_p, grid_u, grid_u_n = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n =
        dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n
    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half,
    dyn_data.grid_lnp_half,
    dyn_data.grid_p_full,
    dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps = dyn_data.grid_lnps

    grid_div, grid_absvor, grid_vor =
        dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half =
        dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half

    grid_energy_full, spe_energy = dyn_data.grid_energy_full, dyn_data.spe_energy

    grid_q_n = dyn_data.grid_q_n
    grid_q_c = dyn_data.grid_q_c
    grid_q_p = dyn_data.grid_q_p

    grid_δq = dyn_data.grid_δq

    grid_virtual_t = vert_coord.virtual_temperature
    grid_w_full = dyn_data.grid_w_full

    K_E = dyn_data.K_E
    # pressure difference
    grid_Δp = dyn_data.grid_Δp
    # temporary variables
    grid_δQ = dyn_data.grid_d_full1

    # incremental quantities
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu,
    dyn_data.grid_δv,
    dyn_data.grid_δps,
    dyn_data.grid_δlnps,
    dyn_data.grid_δt
    integrator = semi_implicit.integrator
    Δt = Get_Δt(integrator)

    # Dynamics conserves the previous time-level invariants.
    mean_ps_p, mean_energy_p, mean_moisture_p = Compute_Corrections_Init(
        mesh,
        vert_coord,
        atmo_data,
        grid_ps_p,
        grid_energy_full,
        grid_u_p,
        grid_v_p,
        grid_t_p,
        grid_q_p,
    )

    # --- Pressure-related --- #
    # Updates pressure variables and computes gradients ∇pₛ.
    # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
    Pressure_Variables!(
        vert_coord,
        grid_ps,
        grid_p_half,
        grid_Δp,
        grid_lnp_half,
        grid_p_full,
        grid_lnp_full,
    )

    # compute ∇ps = ∇lnps * ps
    Compute_Gradients!(mesh, spe_lnps_c, grid_dλ_ps, grid_dθ_ps)
    grid_dλ_ps .*= grid_ps
    grid_dθ_ps .*= grid_ps
    Compute_Virtual_Temperature!(grid_virtual_t, atmo_data, grid_t, grid_q_c)
    # --- Pressure-related --- #

    # --- Adiabatic Tendencies --- #
    # Compute vertical velocity and moist equation-of-state pressure work.
    # Pressure gradient forces are included here.
    Four_In_One!(
        vert_coord,
        atmo_data,
        grid_div,
        grid_u,
        grid_v,
        grid_ps,
        grid_Δp,
        grid_lnp_half,
        grid_lnp_full,
        grid_p_full,
        grid_dλ_ps,
        grid_dθ_ps,
        grid_virtual_t,
        grid_M_half,
        grid_w_full,
        grid_δu,
        grid_δv,
        grid_δps,
        grid_δt,
    )

    Compute_Geopotential!(
        vert_coord,
        atmo_data,
        grid_lnp_half,
        grid_lnp_full,
        grid_virtual_t,
        grid_geopots,
        grid_geopot_full,
        grid_geopot_half,
    )

    grid_δlnps .= grid_δps ./ grid_ps
    Trans_Grid_To_Spherical!(mesh, grid_δlnps, spe_δlnps)
    # --- Adiabatic Tendencies --- #

    # --- Advections --- #
    # Computes non-linear advection terms.
    Vert_Advection!(
        vert_coord,
        grid_u,
        grid_Δp,
        grid_M_half,
        Δt,
        vert_coord.vert_advect_scheme,
        grid_δQ,
    )
    grid_δu .+= grid_δQ
    Vert_Advection!(
        vert_coord,
        grid_v,
        grid_Δp,
        grid_M_half,
        Δt,
        vert_coord.vert_advect_scheme,
        grid_δQ,
    )
    grid_δv .+= grid_δQ
    Vert_Advection!(
        vert_coord,
        grid_t,
        grid_Δp,
        grid_M_half,
        Δt,
        vert_coord.vert_advect_scheme,
        grid_δQ,
    )
    grid_δt .+= grid_δQ

    # Moisture is a grid-only prognostic. The finite-volume transport starts
    # from q_c during startup and q_p during leapfrog, then applies H followed
    # by V over the same effective interval as the dry dynamical core.
    if config.moisture_processes
        grid_q_base = integrator.init_step ? grid_q_c : grid_q_p
        Advance_Grid_Tracer!(
            dyn_data.tracer_workspace,
            mesh,
            grid_q_n,
            grid_q_base,
            grid_u,
            grid_v,
            grid_Δp,
            grid_M_half,
            Δt,
        )
        @. grid_δq = (grid_q_n - grid_q_base) / Δt
        _filter_grid_humidity!(
            integrator,
            mesh,
            vert_coord,
            atmo_data,
            grid_q_p,
            grid_q_c,
            grid_q_n,
            grid_ps,
        )
    else
        copyto!(grid_q_n, integrator.init_step ? grid_q_c : grid_q_p)
    end

    # temperature
    Add_Horizontal_Advection!(mesh, spe_t_c, grid_u, grid_v, grid_δt)
    Trans_Grid_To_Spherical!(mesh, grid_δt, spe_δt)
    # --- Advections --- #

    # --- Vector Invariant Formulation --- #
    # Compute vorticity and divergence tendencies (ηv, ηu) and kinetic energy (E = Φ + ½(u² + v²)).
    grid_absvor = dyn_data.grid_absvor
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)

    grid_δu .+= grid_absvor .* grid_v
    grid_δv .-= grid_absvor .* grid_u

    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)

    grid_energy_full .= grid_geopot_full .+ 0.5 * (grid_u .^ 2 + grid_v .^ 2)
    Trans_Grid_To_Spherical!(mesh, grid_energy_full, spe_energy)
    Apply_Laplacian!(mesh, spe_energy)
    spe_δdiv .-= spe_energy
    # --- Vector Invariant Formulation --- #

    # --- Implicit Correction --- #
    # Inverts the Helmholtz equation in spectral space to stabilize gravity waves.
    Implicit_Correction!(
        semi_implicit,
        vert_coord,
        atmo_data,
        spe_div_c,
        spe_div_p,
        spe_lnps_c,
        spe_lnps_p,
        spe_t_c,
        spe_t_p,
        spe_δdiv,
        spe_δlnps,
        spe_δt,
    )
    # --- Implicit Corrections --- #

    # --- Spectral Damping --- #
    # Applies ∇²ⁿ hyper-diffusion to dampen small-scale noise.
    Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
    Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
    Compute_Spectral_Damping!(integrator, spe_t_c, spe_t_p, spe_δt)
    # --- Spectral Damping --- #

    # --- Time Integration --- #
    Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
    Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
    Filtered_Leapfrog!(integrator, spe_δlnps, spe_lnps_p, spe_lnps_c, spe_lnps_n)
    Filtered_Leapfrog!(integrator, spe_δt, spe_t_p, spe_t_c, spe_t_n)

    Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_n, grid_lnps)
    grid_ps_n .= exp.(grid_lnps)
    Trans_Spherical_To_Grid!(mesh, spe_t_n, grid_t_n)
    # --- Time Integration --- #

    # --- Conservation Fixer --- #
    # Corrects the final predicted state to match global integrals.
    Compute_Corrections!(
        mesh,
        vert_coord,
        atmo_data,
        mean_ps_p,
        grid_ps_n,
        spe_lnps_n,
        mean_energy_p,
        grid_energy_full,
        grid_u_n,
        grid_v_n,
        grid_t_n,
        spe_t_n,
        mean_moisture_p,
        grid_q_n,
    )
    # --- Conservation Fixer --- #

    # Keep the unfiltered grid current as the conservation ledger that rotates
    # into grid_*_p. The Asselin-filtered spectral current is numerical history,
    # not a second authoritative physical state; reconstructing it on the grid
    # makes the global fixer conserve the dissipated computational mode and is
    # unstable for this leapfrog formulation.

    if advance_time
        Time_Advance!(dyn_data)
        Pressure_Variables!(
            vert_coord,
            grid_ps,
            grid_p_half,
            grid_Δp,
            grid_lnp_half,
            grid_p_full,
            grid_lnp_full,
        )
    end

    return
end



function Get_Topography!(grid_geopots::Array{Float64,3})

    grid_geopots .= 0.0

end



"""
    Spectral_Initialize_Fields!(
        config,
        mesh, vert_coord, 
        atmo_data, dyn_data,
        sea_level_ps_ref, init_t, grid_geopots
    )

Initializes the prognostic variables (u, v, T, ln(pₛ)) for the start of the simulation. 
It sets up a balanced initial state (often a zonal jet or a state of rest) perturbed 
by small-amplitude noise in the vorticity field to initiate baroclinic wave growth.

### Parameters
    - config: Model configuration struct.
    - mesh: Spectral spherical mesh properties.
    - vert_coord: Vertical coordinate definitions.
    - atmo_data: Atmospheric constants.
    - dyn_data: Data container for state variables. Modified in-place.
    - sea_level_ps_ref: Reference sea-level pressure (Pa).
    - init_t: Initial isothermal temperature (K) used for the basic state.
    - grid_geopots: Surface geopotential (topography) [nλ, nθ, 1].

### Returns
    - nothing

### Modified
    - dyn_data

### Notes
    - **Surface Pressure Initialization:**
    Derived hydrostatically from the surface geopotential (`grid_geopots`) and the 
    reference sea-level pressure, assuming an isothermal atmosphere at `init_t`:
    ln(pₛ) = ln(pₛ_ref) - Φ_surf / (R_d ⋅ T_init)

    - **Perturbation:**
    Adds a small Gaussian noise or specific spectral mode perturbation to the 
    vorticity field (`spe_vor_c`) to break symmetry and allow baroclinic instability 
    to develop. The perturbation magnitude is scaled by `1.0e-7/sqrt(2.0)`.

    - **Consistency:**
    After setting the spectral vorticity and divergence, it calls `UV_Grid_From_Vor_Div!` 
    and `Pressure_Variables!` to ensure that the grid-point winds (u, v) and 
    diagnostic pressure variables (p_half, Δp, etc.) are fully consistent with 
    the spectral state.
    
    - **Moisture:**
    This legacy initializer is dry. Moist benchmark humidity is initialized by
    `Initialize_Atmos_State!`.

"""
function Spectral_Initialize_Fields!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    sea_level_ps_ref::Float64,
    init_t::Float64,
    grid_geopots::Array{Float64,3},
)

    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
    rdgas = atmo_data.rdgas

    spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c =
        dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
    spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p =
        dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p

    grid_u, grid_v, grid_ps, grid_t =
        dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
    grid_u_p, grid_v_p, grid_ps_p, grid_t_p =
        dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p

    grid_lnps, grid_vor, grid_div = dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div

    grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half,
    dyn_data.grid_Δp,
    dyn_data.grid_lnp_half,
    dyn_data.grid_p_full,
    dyn_data.grid_lnp_full

    grid_t .= init_t

    # dΦ/dlnp = -RT    Δp = -ΔΦ/RT
    grid_lnps[:, :, 1] .= log(sea_level_ps_ref) .- grid_geopots[:, :, 1] ./ (rdgas * init_t)
    grid_ps .= exp.(grid_lnps)


    spe_div_c .= 0.0
    spe_vor_c .= 0.0

    # # initial perturbation
    num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical

    initial_perturbation = 1.0e-7 / sqrt(2.0)
    # initial vorticity perturbation used in benchmark code
    # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
    for k = nd-2:nd
        spe_vor_c[2, 5, k] = initial_perturbation
        spe_vor_c[6, 9, k] = initial_perturbation
        spe_vor_c[2, 4, k] = initial_perturbation
        spe_vor_c[6, 8, k] = initial_perturbation
    end

    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)


    # initial spectral fields (and spectrally-filtered) grid fields

    Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
    Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)

    Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps_c)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_c, grid_lnps)
    grid_ps .= exp.(grid_lnps)


    Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)

    UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)

    Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)

    #update pressure variables for hs forcing
    Pressure_Variables!(
        vert_coord,
        grid_ps,
        grid_p_half,
        grid_Δp,
        grid_lnp_half,
        grid_p_full,
        grid_lnp_full,
    )


    spe_vor_p .= spe_vor_c
    spe_div_p .= spe_div_c
    spe_lnps_p .= spe_lnps_c
    spe_t_p .= spe_t_c


    grid_u_p .= grid_u
    grid_v_p .= grid_v
    grid_ps_p .= grid_ps
    grid_t_p .= grid_t

    # This legacy dry initializer is retained only for the dry public API. The
    # moist benchmark is initialized by `Initialize_Atmos_State!`.
    fill!(dyn_data.grid_q_c, 0.0)
    copyto!(dyn_data.grid_q_p, dyn_data.grid_q_c)

end



function Atmosphere_Update!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    semi_implicit::Semi_Implicit_Solver,
    physics_params::Dict{String,Any},
)
    # First form a corrected dynamics-only provisional next state. Physics is
    # then a direct time-split map of that state, never a leapfrog tendency.
    Spectral_Dynamics!(
        config,
        mesh,
        vert_coord,
        atmo_data,
        dyn_data,
        semi_implicit,
        advance_time = false,
    )

    Pressure_Variables!(
        vert_coord,
        dyn_data.grid_ps_n,
        dyn_data.grid_p_half,
        dyn_data.grid_Δp,
        dyn_data.grid_lnp_half,
        dyn_data.grid_p_full,
        dyn_data.grid_lnp_full,
    )
    Spectral_Physics!(
        config,
        mesh,
        vert_coord,
        atmo_data,
        dyn_data,
        semi_implicit,
        physics_params,
    )
    physics_changed_state = any(
        get(physics_params, key, false) for key in (
            "do_Betts_Miller",
            "do_Lscale_Cond",
            "do_Sensible_Heating",
            "do_Surface_Evaporation",
            "do_Implicit_PBL_Scheme",
            "do_HS_Forcing",
            "do_LRF",
        )
    )
    physics_changed_state &&
        _synchronize_physics_next!(mesh, vert_coord, atmo_data, dyn_data)
    Time_Advance!(dyn_data)

    grid_ps, grid_Δp, grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full =
        dyn_data.grid_ps_c,
        dyn_data.grid_Δp,
        dyn_data.grid_p_half,
        dyn_data.grid_lnp_half,
        dyn_data.grid_p_full,
        dyn_data.grid_lnp_full

    grid_t = dyn_data.grid_t_c
    grid_geopots, grid_z_full, grid_z_half =
        dyn_data.grid_geopots, dyn_data.grid_z_full, dyn_data.grid_z_half

    grid_q_c = dyn_data.grid_q_c

    Compute_Pressures_And_Heights!(
        atmo_data,
        vert_coord,
        grid_ps,
        grid_geopots,
        grid_t,
        grid_p_half,
        grid_Δp,
        grid_lnp_half,
        grid_p_full,
        grid_lnp_full,
        grid_z_full,
        grid_z_half,
        grid_q_c,
    )

    return
end

end
