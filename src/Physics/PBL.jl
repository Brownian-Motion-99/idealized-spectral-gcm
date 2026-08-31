using Base.Threads
using ...Atmo_Data_Module
using ...Dyn_Data_Module
using ...Spectral_Spherical_Mesh_Module

struct PBL_Workspace
    V_c::Matrix{Float64}
    za::Matrix{Float64}
end

PBL_Workspace(nλ::Int, nθ::Int) =
    PBL_Workspace(zeros(nλ, nθ), zeros(nλ, nθ))


"""
    Calculate_V_c_za!(
        workspace,
        atmo_data,
        grid_p_half, grid_ps,
        grid_u, grid_v,
        grid_t, grid_q
    )

Calculates the surface wind speed (`V_c`) and hydrostatic height of the lowest
model level (`za`) needed by the MITC bulk tendencies.

### Parameters
    - atmo_data: Structure containing physical constants (R_d, g).
    
    - grid_p_half: Pressure at layer interfaces [nλ, nθ, nd+1].
    - grid_ps: Surface pressure [nλ, nθ, 1].
    
    - grid_u, grid_v: Zonal and meridional wind components [nλ, nθ, nd].
    
    - grid_t: Temperature [nλ, nθ, nd].
    - grid_q: Specific humidity [nλ, nθ, nd].

### Returns
    - V_c: Magnitude of the horizontal wind at the lowest model level [nλ, nθ].
    - za: Geometric height of the lowest model level (anemometer height proxy) [nλ, nθ],
      calculated from the surface pressure and the interface above the lowest layer
      following Thatcher and Jablonowski (2016), Eq. (12).
"""
function Calculate_V_c_za!(
    workspace::PBL_Workspace,
    atmo_data::Atmo_Data,
    grid_p_half::Array{Float64,3},
    grid_ps::Array{Float64,3},
    grid_u::Array{Float64,3},
    grid_v::Array{Float64,3},
    grid_t::Array{Float64,3},
    grid_q::AbstractArray{Float64,3},
)

    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64
    Rd, Rv, grav =
        atmo_data.rdgas::Float64, atmo_data.rvgas::Float64, atmo_data.grav::Float64

    Rd_g = Rd / grav
    virtual_coefficient = Rv / Rd - 1.0
    V_c, za = workspace.V_c, workspace.za

    @threads for j = 1:nθ
        for i = 1:nλ

            # Unpack
            us_val = grid_u[i, j, nd]
            vs_val = grid_v[i, j, nd]
            ts_val = grid_t[i, j, nd]
            qs_val = grid_q[i, j, nd]
            p_above_val = grid_p_half[i, j, nd]
            ps_val = grid_ps[i, j, 1]

            # tv
            tv_factor = 1.0 + virtual_coefficient * qs_val

            # Surface wind speed
            V_c[i, j] = sqrt(us_val^2 + vs_val^2)

            # Height of the lowest full level. Equation (12) of Thatcher and
            # Jablonowski (2016) uses the interface between the two lowest full
            # levels; the factor 1/2 places the full level halfway through the
            # lowest layer in log-pressure height.
            tvs = ts_val * tv_factor
            za[i, j] = Rd_g * tvs * log(ps_val / p_above_val) * 0.5

        end
    end

    return V_c, za
end



"""
    Sensible_Heating!(
        mesh, atmo_data,
        grid_p_half, grid_t, grid_shflx,
        V_c, za,
        Δt, 
        C_H,
        lower_boundary_temperature
    )

Computes the exchange of sensible heat between the surface and the lowest atmospheric 
layer. The implementation uses a backward-implicit update to solve for the new 
temperature, ensuring stability without restricting the model time step.

### Parameters
    - mesh: Spectral mesh properties (provides longitude λc and latitude θc).
    - atmo_data: Atmospheric constants (cₚ, nλ, nθ).
    - grid_p_half: Layer-interface pressures used to obtain the exact mass
      `Δp_bottom/g` of the lowest model layer.
    
    - grid_t: Temperature field [nλ, nθ, nd]. Modified in-place (lowest level updated).
    - grid_shflx: Output sensible heat flux [nλ, nθ, 1] in W/m². Modified in-place.
    
    - V_c: Surface wind speed magnitude [nλ, nθ].
    - za: Height of the lowest model level [nλ, nθ].
    
    - Δt: Physics time step.
    
    - C_H: Bulk aerodynamic transfer coefficient for heat (Stanton number). Default 0.0044.
    - lower_boundary_temperature: Function `(longitude, latitude) -> temperature_K`.
      Defaults to `Default_Lower_Boundary_Temperature`.

### Returns
    - nothing

### Modified
    - grid_t (at index k=nd)
    - grid_shflx

"""
function Sensible_Heating!(
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
    grid_p_half::Array{Float64,3},
    grid_t::Array{Float64,3},
    grid_shflx::Array{Float64,3},
    V_c::Array{Float64,2},
    za::Array{Float64,2},
    Δt::Int64,
    C_H::Float64 = 0.0044,
    lower_boundary_temperature = Default_Lower_Boundary_Temperature,
)

    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64
    cp = atmo_data.cp_air::Float64
    grav = atmo_data.grav::Float64
    λc, θc = mesh.λc, mesh.θc
    size(grid_p_half) == (nλ, nθ, nd + 1) ||
        throw(DimensionMismatch("sensible-heating p_half has incorrect size"))

    @threads for j = 1:nθ

        for i = 1:nλ
            surface_temperature = Float64(lower_boundary_temperature(λc[i], θc[j]))

            # Implicit coef.
            # λ = (C_H * |V| * Δt) / za
            lambda = (C_H * V_c[i, j] * Float64(Δt)) / za[i, j]

            # New temperature
            # T_new = (T_old + λ * T_s) / (1 + λ)
            t_old = grid_t[i, j, nd]
            t_new = (t_old + lambda * surface_temperature) / (1 + lambda)
            grid_t[i, j, nd] = t_new

            # Diagnose the exact finite-volume atmospheric energy increment.
            # The MITC tendency still uses za; only its reported flux uses the
            # actual pressure-coordinate mass of the bottom model layer.
            Δp_bottom = grid_p_half[i, j, nd+1] - grid_p_half[i, j, nd]
            isfinite(Δp_bottom) && Δp_bottom > 0.0 ||
                throw(DomainError(Δp_bottom, "bottom-layer pressure thickness must be positive"))
            grid_shflx[i, j, 1] =
                (Δp_bottom / grav) * cp * (t_new - t_old) / Float64(Δt)

        end

    end

end



"""
    Surface_Evaporation!(
        mesh, atmo_data,
        grid_ps, grid_p_half,
        grid_q, grid_lhflx,
        V_c, za,
        Δt, 
        C_E,
        lower_boundary_temperature
    )

Computes the vertical transport of water vapor from the surface to the lowest atmospheric 
layer. The scheme employs a backward-implicit time integration to ensure numerical stability 
and enforces a "no-dew" condition (evaporation only).

### Parameters
    - mesh: Spectral mesh properties (provides longitude λc and latitude θc).
    - atmo_data: Atmospheric constants (Lᵥ, Rᵥ, etc.).
    
    - grid_ps: Surface pressure [nλ, nθ, 1].
    - grid_p_half: Layer-interface pressures used to obtain the exact mass
      `Δp_bottom/g` of the lowest model layer.
    
    - grid_q: Specific humidity field [nλ, nθ, nd]. Modified in-place (lowest level updated).
    - grid_lhflx: Output latent heat flux [nλ, nθ, 1] in W/m². Modified in-place.
    
    - V_c: Surface wind speed magnitude [nλ, nθ].
    - za: Height of the lowest model level [nλ, nθ].
    
    - Δt: Physics time step.
    
    - C_E: Bulk aerodynamic transfer coefficient for moisture (Dalton number). Default 0.0044.
    - lower_boundary_temperature: Function `(longitude, latitude) -> temperature_K`.
      Defaults to `Default_Lower_Boundary_Temperature`.

### Returns
    - nothing

### Modified
    - grid_q (at index k=nd)
    - grid_lhflx

"""
function Surface_Evaporation!(
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
    grid_ps::Array{Float64,3},
    grid_p_half::Array{Float64,3},
    grid_q::AbstractArray{Float64,3},
    grid_lhflx::Array{Float64,3},
    V_c::Array{Float64,2},
    za::Array{Float64,2},
    Δt::Int64,
    C_E::Float64 = 0.0044,
    lower_boundary_temperature = Default_Lower_Boundary_Temperature,
)

    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64
    Lv = atmo_data.Lv::Float64
    grav = atmo_data.grav::Float64
    epsilon = atmo_data.rdgas / atmo_data.rvgas
    λc, θc = mesh.λc, mesh.θc
    size(grid_p_half) == (nλ, nθ, nd + 1) ||
        throw(DimensionMismatch("surface-evaporation p_half has incorrect size"))

    @threads for j = 1:nθ

        for i = 1:nλ
            surface_temperature = Float64(lower_boundary_temperature(λc[i], θc[j]))

            # Unpack
            ps_val = grid_ps[i, j, 1]
            qs_val = grid_q[i, j, nd]

            # Surface saturated specific humidity
            qs = Saturation_Specific_Humidity(surface_temperature, ps_val, epsilon)

            # Implicit coef.
            # λ = (C_E * |V| * Δt) / za
            lambda = (C_E * V_c[i, j] * Float64(Δt)) / za[i, j]

            # New specific humidity
            # q_new = (q_old + λ * q_target) / (1 + λ)
            q_target = max(qs_val, qs)
            q_old = qs_val
            q_new = (q_old + lambda * q_target) / (1 + lambda)
            grid_q[i, j, nd] = q_new

            # Diagnose the exact finite-volume atmospheric water increment in
            # latent-energy units. Dry_Air_Adjustment! subsequently preserves
            # this pre-adjustment layer-water amount exactly.
            Δp_bottom = grid_p_half[i, j, nd+1] - grid_p_half[i, j, nd]
            isfinite(Δp_bottom) && Δp_bottom > 0.0 ||
                throw(DomainError(Δp_bottom, "bottom-layer pressure thickness must be positive"))
            grid_lhflx[i, j, 1] =
                (q_new - q_old) * (Δp_bottom / grav) * Lv / Float64(Δt)

        end

    end

end



abstract type PBLTop end
struct PressureLevelBasedPBLTop <: PBLTop end
struct ModelLevelBasedPBLTop <: PBLTop end

function SelectPBLTop(sym::Symbol, value::Union{Int64,Float64})
    if sym == :PressureLevel
        @assert isa(value, Float64) "pbl_top_mode $sym must have Float64 pbl_top_value given."
        return PressureLevelBasedPBLTop(), value

    elseif sym == :ModelLevel
        @assert isa(value, Int64) "pbl_top_mode $sym must have Int64 pbl_top_value given."
        return ModelLevelBasedPBLTop(), value

    else
        error("Implicit PBL Mixing Scheme: Unknown PBLTop Symbol: $sym")
    end
end

PBL_Top_Symbol(::PressureLevelBasedPBLTop) = :PressureLevel
PBL_Top_Symbol(::ModelLevelBasedPBLTop) = :ModelLevel



"""
    Implicit_PBL_Mixing!(
        atmo_data,
        grid_p_full, grid_p_half, 
        grid_t, grid_q,
        K_E, 
        V_c, za, rho,
        physics_params,
        Δt, 
        C_D
    )

### Parameters
    - atmo_data: Atmospheric constants (R, cₚ, g).
    
    - grid_p_full: Pressure at layer centers.
    - grid_p_half: Pressure at layer interfaces.
    
    - grid_t: Temperature field [nλ, nθ, nd]. Modified in-place.
    - grid_q: Specific humidity field [nλ, nθ, nd]. Modified in-place.
    
    - K_E: Turbulent diffusivity array [nλ, nθ, nd+1]. Modified in-place (diagnostic).
    
    - V_c: Surface wind speed magnitude.
    - za: Height of the lowest model level.
    - physics_params: Dictionary defining the PBL top definition (:PressureLevel or :ModelLevel).
    
    - Δt: Physics time step.
    
    - C_D: Drag coefficient for momentum/heat. Default 0.0044.

### Returns
    - nothing

### Modified
    - grid_t
    - grid_q
    - K_E

This function includes:
    1. Calculating turbulent diffusivity (K_E)
    2. Implicitly calculating the temperature and moisture at the next timestep

Implicit Euler scheme for tracer concentration (ϕ), trancers can be moisture, potential temperature, etc.:
Δϕ^{n}_{k} = (ϕ^{n+1}_{k} - ϕ^{n}_k) / Δt

Rewrite in flux form: 
Δϕ^{n}_{k} = (Flux^{n+1}_{k+1/2} - Flux^{n+1}_{k-1/2}) / Δp_{k}

The flux is given by:
Flux^{n+1}_{k+1/2} ≈ ((g^2 * ρ^2 * K_E) / Δp_{k+1/2}) * (ϕ^{n+1}_{k+1} - ϕ^{n+1}_{k})

Replace the flux with three ϕs at the next timestep:
ϕ^{n}_{k} = -CA * ϕ^{n+1}_{k-1} + (1 + CA + CC) * ϕ^{n+1}_{k} - CC * ϕ^{n+1}_{k+1}

This equation is solved using Thomas Algorithm

see also: Reed and Jablonowski (JAMES, 2012)

"""
function Implicit_PBL_Mixing!(
    atmo_data::Atmo_Data,
    grid_p_full::Array{Float64,3},
    grid_p_half::Array{Float64,3},
    grid_t::Array{Float64,3},
    grid_q::Array{Float64,3},
    K_E::Array{Float64,3},
    V_c::Array{Float64,2},
    za::Array{Float64,2},
    physics_params::Dict{String,Any},
    Δt::Int64,
    C_D::Float64 = 0.0044,
)

    pbl_top_mode, pbl_top_value = SelectPBLTop(
        physics_params["PBL_Top_Mode"]::Symbol,
        physics_params["PBL_Top_Value"]::Union{Int64,Float64},
    )

    nλ, nθ, nd = atmo_data.nλ, atmo_data.nθ, atmo_data.nd
    Rd, cp, grav = atmo_data.rdgas, atmo_data.cp_air, atmo_data.grav
    virtual_coefficient = atmo_data.rvgas / Rd - 1.0

    p_scale = 10000.0
    p0 = 100000.0
    Rd_cp = Rd / cp
    grav_sq = grav^2

    @threads for j = 1:nθ

        # --- Local buffers for threads --- #
        CA = Vector{Float64}(undef, nd)      # Coupling coef. to the layer above
        CC = Vector{Float64}(undef, nd)      # Coupling coef. to the layer below
        CE = Vector{Float64}(undef, nd + 1)    # Eliminator
        CF = Vector{Float64}(undef, nd + 1)    # Source of latent heat
        CFt = Vector{Float64}(undef, nd + 1)    # Source of sensible heat
        # --- Local buffers for threads --- #

        for i = 1:nλ

            # --- Load data --- #
            V_c_val = V_c[i, j]
            za_val = za[i, j]
            # --- Load data --- #

            # --- Calculate K_E --- #
            # Within PBL: K_E = C_D * Vs * za
            # Above PBL:  K_E = C_D * Vs * za * exp(-((p_pbl_top - p) / H)^2)
            if isa(pbl_top_mode, PressureLevelBasedPBLTop)
                p_pbl_top = pbl_top_value
                for k = 1:nd+1
                    if grid_p_half[i, j, k] >= p_pbl_top
                        K_E[i, j, k] = C_D * V_c_val * za_val
                    else
                        K_E[i, j, k] =
                            C_D *
                            V_c_val *
                            za_val *
                            exp(-((p_pbl_top - grid_p_half[i, j, k]) / p_scale)^2)
                    end
                end

            elseif isa(pbl_top_mode, ModelLevelBasedPBLTop)
                p_end_constant_mixing = 85000.0
                idk_pbl_top = pbl_top_value
                for k = nd+1-idk_pbl_top:nd+1
                    K_E[i, j, k] = C_D * V_c_val * za_val
                end
                for k = 1:nd+1-idk_pbl_top-1
                    K_E[i, j, k] =
                        C_D *
                        V_c_val *
                        za_val *
                        exp(-((p_end_constant_mixing - grid_p_half[i, j, k]) / p_scale)^2)
                end

            else
                error("Implicit PBL Mixing Scheme: Unknown PBLTop Symbol: $pbl_top_mode")
            end
            # --- Calculate K_E --- #

            # --- Finite volume coupling coef. --- #
            # CA = Δt * (g^2 * ρ^2 * K_E) / (Δp_{k+1/2} * Δp_{k})
            # CC = Δt * (g^2 * ρ^2 * K_E) / (Δp_{k-1/2} * Δp_{k})
            CA[nd] = 0.0
            CC[1] = 0.0
            for k = 1:nd-1
                rpdel_k = 1 / (grid_p_half[i, j, k+1] - grid_p_half[i, j, k])
                rpdel_kp1 = 1 / (grid_p_half[i, j, k+2] - grid_p_half[i, j, k+1])
                tv_upper = grid_t[i, j, k] *
                           (1.0 + virtual_coefficient * grid_q[i, j, k])
                tv_lower = grid_t[i, j, k+1] *
                           (1.0 + virtual_coefficient * grid_q[i, j, k+1])
                rho_interface =
                    grid_p_half[i, j, k+1] / (Rd * 0.5 * (tv_upper + tv_lower))
                isfinite(rho_interface) && rho_interface > 0 || throw(
                    ArgumentError("PBL interface density must be positive and finite"),
                )
                CA[k] =
                    rpdel_k * Float64(Δt) * grav_sq * K_E[i, j, k+1] * rho_interface^2 /
                    (grid_p_full[i, j, k+1] - grid_p_full[i, j, k])
                CC[k+1] =
                    rpdel_kp1 * Float64(Δt) * grav_sq * K_E[i, j, k+1] * rho_interface^2 /
                    (grid_p_full[i, j, k+1] - grid_p_full[i, j, k])
            end
            # --- Finite volume coupling coef. --- #

            # --- Forward sweep --- #
            # Thomas Algorithm
            CE[1] = 0.0
            CE[nd+1] = 0.0
            CF[nd+1] = 0.0
            CFt[nd+1] = 0.0
            for k = nd:-1:1
                CE[k] = CC[k] / (1.0 + CA[k] + CC[k] - CA[k] * CE[k+1])
                CF[k] = (
                    (grid_q[i, j, k] + CA[k] * CF[k+1]) /
                    (1.0 + CA[k] + CC[k] - CA[k] * CE[k+1])
                )
                CFt[k] = (
                    (
                        (p0 / grid_p_full[i, j, k])^(Rd_cp) * grid_t[i, j, k] +
                        CA[k] * CFt[k+1]
                    ) / (1.0 + CA[k] + CC[k] - CA[k] .* CE[k+1])
                )
            end
            # --- Forward sweep --- #

            # --- Backward substitute --- #
            # Loop downward to calculate the moisture and temperature at the next timestep
            # Note that the tracer for temperature is potential temperature,
            # so we have to transform the potential temperature to temperature with Poisson's equation
            grid_q[i, j, 1] = CF[1]
            grid_t[i, j, 1] = CFt[1] * (grid_p_full[i, j, 1] / p0)^Rd_cp
            for k = 2:nd
                grid_q[i, j, k] = CE[k] * grid_q[i, j, k-1] + CF[k]
                grid_t[i, j, k] =
                    (
                        CE[k] * grid_t[i, j, k-1] * (p0 / grid_p_full[i, j, k-1])^Rd_cp +
                        CFt[k]
                    ) * (grid_p_full[i, j, k] / p0)^Rd_cp
            end
            # --- Backward substitute --- #

        end
    end

end
