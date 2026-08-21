using Base.Threads
using ...Atmo_Data_Module
using ...Dyn_Data_Module
using ...Spectral_Spherical_Mesh_Module



"""
    Calculate_V_c_za_rho(
        atmo_data,
        grid_p_half, grid_p_full, grid_ps,
        grid_u, grid_v,
        grid_t, grid_q
    )

Calculates the surface wind speed (V_c), the height of the lowest model layer (za), 
and the air density profile (ρ) needed for bulk aerodynamic flux calculations and 
vertical diffusion.

### Parameters
    - atmo_data: Structure containing physical constants (R_d, g).
    
    - grid_p_half: Pressure at layer interfaces [nλ, nθ, nd+1].
    - grid_p_full: Pressure at layer centers [nλ, nθ, nd].
    - grid_ps: Surface pressure [nλ, nθ, 1].
    
    - grid_u, grid_v: Zonal and meridional wind components [nλ, nθ, nd].
    
    - grid_t: Temperature [nλ, nθ, nd].
    - grid_q: Specific humidity [nλ, nθ, nd].

### Returns
    - V_c: Magnitude of the horizontal wind at the lowest model level [nλ, nθ].
    - za: Geometric height of the lowest model level (anemometer height proxy) [nλ, nθ],
      calculated from the surface pressure and the interface above the lowest layer
      following Thatcher and Jablonowski (2016), Eq. (12).
    - rho: Air density at layer centers [nλ, nθ, nd].

"""
function Calculate_V_c_za_rho(
    atmo_data::Atmo_Data,
    grid_p_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_ps::Array{Float64, 3},
    grid_u::Array{Float64, 3}, grid_v::Array{Float64, 3},
    grid_t::Array{Float64, 3}, grid_q::AbstractArray{Float64, 3}
)
    
    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64
    Rd, grav   = atmo_data.rdgas::Float64, atmo_data.grav::Float64

    Rd_g   = Rd / grav
    inv_Rd = 1.0 / Rd

    # Allocation
    V_c = zeros(Float64, nλ, nθ)
    za  = zeros(Float64, nλ, nθ)
    rho = zeros(Float64, nλ, nθ, nd)

    @threads for j = 1:nθ
        for i = 1:nλ

            # Unpack
            us_val      = grid_u[i, j, nd]
            vs_val      = grid_v[i, j, nd]
            ts_val      = grid_t[i, j, nd]
            qs_val      = grid_q[i, j, nd]
            p_above_val = grid_p_half[i, j, nd]
            ps_val      = grid_ps[i, j, 1]
            
            # tv
            tv_factor = 1.0 + 0.608 * qs_val

            # Surface wind speed
            V_c[i, j] = sqrt(us_val^2 + vs_val^2)

            # Height of the lowest full level. Equation (12) of Thatcher and
            # Jablonowski (2016) uses the interface between the two lowest full
            # levels; the factor 1/2 places the full level halfway through the
            # lowest layer in log-pressure height.
            tvs      = ts_val * tv_factor
            za[i, j] = Rd_g * tvs * log(ps_val / p_above_val) * 0.5

            # Density uses the virtual temperature of each level, not the
            # lowest-level humidity for the entire column.
            for k = 1:nd
                p_full_val = grid_p_full[i, j, k]
                t_val      = grid_t[i, j, k]
                q_val      = grid_q[i, j, k]
                tv_factor_k = 1.0 + 0.608 * q_val
                rho[i, j, k] = p_full_val * inv_Rd / (t_val * tv_factor_k)
            end
            
        end
    end
    
    return V_c, za, rho
end



"""
    Sensible_Heating!(
        mesh, atmo_data,
        grid_t, grid_shflx,
        V_c, za, rho,
        Δt, 
        C_H
    )

Computes the exchange of sensible heat between the surface and the lowest atmospheric 
layer. The implementation uses a backward-implicit update to solve for the new 
temperature, ensuring stability without restricting the model time step.

### Parameters
    - mesh: Spectral mesh properties (provides latitude θc for Held-Suarez forcing).
    - atmo_data: Atmospheric constants (cₚ, nλ, nθ).
    
    - grid_t: Temperature field [nλ, nθ, nd]. Modified in-place (lowest level updated).
    - grid_shflx: Output sensible heat flux [nλ, nθ, 1] in W/m². Modified in-place.
    
    - V_c: Surface wind speed magnitude [nλ, nθ].
    - za: Height of the lowest model level [nλ, nθ].
    - rho: Air density at the lowest model level [nλ, nθ, nd].
    
    - Δt: Physics time step.
    
    - C_H: Bulk aerodynamic transfer coefficient for heat (Stanton number). Default 0.0044.

### Returns
    - nothing

### Modified
    - grid_t (at index k=nd)
    - grid_shflx

"""
function Sensible_Heating!(
    mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_t::Array{Float64, 3}, grid_shflx::Array{Float64, 3},
    V_c::Array{Float64, 2}, za::Array{Float64, 2}, rho::Array{Float64, 3},
    Δt::Int64,
    C_H::Float64=0.0044
)
    
    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64
    cp         = atmo_data.cp_air::Float64
    θc         = mesh.θc

    deg_factor = 26.0 * pi / 180.0
    denom_lat  = 2.0 * deg_factor^2

    @threads for j = 1:nθ

        # Surface temperature (Held-Suarez)
        lat_factor = θc[j]^2 / denom_lat
        tsj        = 29.0 * exp(-lat_factor) + 271.0

        for i = 1:nλ
            
            # Implicit coef.
            # λ = (C_H * |V| * Δt) / za
            lambda = (C_H * V_c[i, j] * Float64(Δt)) / za[i, j]

            # New temperature
            # T_new = (T_old + λ * T_s) / (1 + λ)
            t_old = grid_t[i, j, nd]
            t_new = (t_old + lambda * tsj) / (1 + lambda)
            grid_t[i, j, nd] = t_new

            # Back-calculate sensible heat flux (W/m^2)
            # energy diff. = mass of Col. * cp * dT
            # mass of Col. ≈ rho * za
            # Flux = energy diff. / Δt ≈ rho * za * cp * dT / Δt
            grid_shflx[i, j, 1] = rho[i, j, nd] * za[i, j] * cp * (t_new - t_old) / Float64(Δt)

        end
        
    end

end



"""
    Surface_Evaporation!(
        mesh, atmo_data,
        grid_ps,
        grid_q, grid_lhflx,
        V_c, za, rho,
        Δt, 
        C_E
    )

Computes the vertical transport of water vapor from the surface to the lowest atmospheric 
layer. The scheme employs a backward-implicit time integration to ensure numerical stability 
and enforces a "no-dew" condition (evaporation only).

### Parameters
    - mesh: Spectral mesh properties (provides latitude θc).
    - atmo_data: Atmospheric constants (Lᵥ, Rᵥ, etc.).
    
    - grid_ps: Surface pressure [nλ, nθ, 1].
    
    - grid_q: Specific humidity field [nλ, nθ, nd]. Modified in-place (lowest level updated).
    - grid_lhflx: Output latent heat flux [nλ, nθ, 1] in W/m². Modified in-place.
    
    - V_c: Surface wind speed magnitude [nλ, nθ].
    - za: Height of the lowest model level [nλ, nθ].
    - rho: Air density at the lowest model level [nλ, nθ, nd].
    
    - Δt: Physics time step.
    
    - C_E: Bulk aerodynamic transfer coefficient for moisture (Dalton number). Default 0.0044.

### Returns
    - nothing

### Modified
    - grid_q (at index k=nd)
    - grid_lhflx

"""
function Surface_Evaporation!(
    mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_ps::Array{Float64, 3},
    grid_q::AbstractArray{Float64,3}, grid_lhflx::Array{Float64, 3},
    V_c::Array{Float64, 2}, za::Array{Float64, 2}, rho::Array{Float64, 3},
    Δt::Int64,
    C_E::Float64=0.0044
)
    
    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64
    Lv         = atmo_data.Lv::Float64
    Rv         = atmo_data.rvgas::Float64
    θc         = mesh.θc

    deg_factor = 26.0 * pi / 180.0
    denom_lat  = 2.0 * deg_factor^2
    const_es   = 611.12
    const_q1   = 0.622
    const_q2   = 0.378
    Lv_Rv      = Lv / Rv
    inv_273    = 1.0 / 273.15

    @threads for j = 1:nθ

        # Surface saturated specific humidity
        lat_factor = θc[j]^2 / denom_lat
        tsj        = 29.0 * exp(-lat_factor) + 271.0
        es         = const_es * exp(Lv_Rv * (inv_273 - 1.0 / tsj))

        for i = 1:nλ

            # Unpack
            ps_val = grid_ps[i, j, 1]
            qs_val = grid_q[i, j, nd]
            
            # Surface saturated specific humidity
            qs = (const_q1 * es) / (ps_val - const_q2 * es)

            # Implicit coef.
            # λ = (C_E * |V| * Δt) / za
            lambda = (C_E * V_c[i, j] * Float64(Δt)) / za[i, j]

            # New specific humidity
            # q_new = (q_old + λ * q_target) / (1 + λ)
            q_target = max(qs_val, qs)
            q_old    = qs_val
            q_new    = (q_old + lambda * q_target) / (1 + lambda)
            grid_q[i, j, nd] = q_new

            # Back-calculate latent heat flux (W/m^2)
            # water mass diff. per unit area = q diff. * mass of col.
            # mass of col. ≈ rho * za
            # Flux = q diff. * mass of col. * Lv / Δt ≈ q diff. * rho * za * Lv / Δt
            grid_lhflx[i, j, 1] = (q_new - q_old) * rho[i, j, nd] * za[i, j] * Lv / Float64(Δt)

        end

    end
    
end 



abstract type PBLTop end
struct PressureLevelBasedPBLTop <: PBLTop end
struct ModelLevelBasedPBLTop <: PBLTop end

function SelectPBLTop(sym::Symbol, value::Union{Int64, Float64})
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
PBL_Top_Symbol(::ModelLevelBasedPBLTop)    = :ModelLevel



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
    - rho: Air density profile.
    
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
    grid_p_full::Array{Float64, 3}, grid_p_half::Array{Float64, 3}, 
    grid_t::Array{Float64, 3}, grid_q::Array{Float64, 3},
    K_E::Array{Float64, 3}, 
    V_c::Array{Float64, 2}, za::Array{Float64, 2}, rho::Array{Float64, 3},
    physics_params::Dict{String, Any},
    Δt::Int64, 
    C_D::Float64=0.0044
)
    
    pbl_top_mode, pbl_top_value = SelectPBLTop(
        physics_params["PBL_Top_Mode"]::Symbol, physics_params["PBL_Top_Value"]::Union{Int64, Float64}
    )

    nλ, nθ, nd   = atmo_data.nλ, atmo_data.nθ, atmo_data.nd
    Rd, cp, grav = atmo_data.rdgas, atmo_data.cp_air, atmo_data.grav

    p_scale   = 10000.0
    p0        = 100000.0
    Rd_cp     = Rd / cp
    grav_sq   = grav^2

    @threads for j = 1:nθ
        
        # --- Local buffers for threads --- #
        CA  = Vector{Float64}(undef, nd)      # Coupling coef. to the layer above
        CC  = Vector{Float64}(undef, nd)      # Coupling coef. to the layer below
        CE  = Vector{Float64}(undef, nd+1)    # Eliminator
        CF  = Vector{Float64}(undef, nd+1)    # Source of latent heat
        CFt = Vector{Float64}(undef, nd+1)    # Source of sensible heat
        # --- Local buffers for threads --- #

        for i = 1:nλ

            # --- Load data --- #
            V_c_val = V_c[i, j]
            za_val  = za[i, j]
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
                        K_E[i, j, k] = C_D * V_c_val * za_val * exp(-((p_pbl_top - grid_p_half[i, j, k]) / p_scale)^2)
                    end
                end

            elseif isa(pbl_top_mode, ModelLevelBasedPBLTop)
                p_end_constant_mixing = 85000.0
                idk_pbl_top           = pbl_top_value
                for k = nd+1-idk_pbl_top:nd+1
                    K_E[i, j, k] = C_D * V_c_val * za_val
                end
                for k = 1:nd+1-idk_pbl_top-1
                    K_E[i, j, k] = C_D * V_c_val * za_val * exp(-((p_end_constant_mixing - grid_p_half[i, j, k]) / p_scale)^2)
                end

            else
                error("Implicit PBL Mixing Scheme: Unknown PBLTop Symbol: $pbl_top_mode")
            end
            # --- Calculate K_E --- #

            # --- Finite volume coupling coef. --- #
            # CA = Δt * (g^2 * ρ^2 * K_E) / (Δp_{k+1/2} * Δp_{k})
            # CC = Δt * (g^2 * ρ^2 * K_E) / (Δp_{k-1/2} * Δp_{k})
            CA[nd] = 0.
            CC[1]  = 0.
            for k = 1:nd-1
                rpdel_k   = 1 / (grid_p_half[i, j, k+1] - grid_p_half[i, j, k])
                rpdel_kp1 = 1 / (grid_p_half[i, j, k+2] - grid_p_half[i, j, k+1])
                CA[k]     = rpdel_k   * Float64(Δt) * grav_sq * K_E[i, j, k+1] * rho[i, j, k+1]^2 / (grid_p_full[i, j, k+1] - grid_p_full[i, j, k])
                CC[k+1]   = rpdel_kp1 * Float64(Δt) * grav_sq * K_E[i, j, k+1] * rho[i, j, k+1]^2 / (grid_p_full[i, j, k+1] - grid_p_full[i, j, k])
            end
            # --- Finite volume coupling coef. --- #

            # --- Forward sweep --- #
            # Thomas Algorithm
            CE[1]     = 0.
            CE[nd+1]  = 0.
            CF[nd+1]  = 0.
            CFt[nd+1] = 0.
            for k = nd:-1:1
                CE[k]  = CC[k] / (1. + CA[k] + CC[k] - CA[k] * CE[k+1])
                CF[k]  = ((grid_q[i, j, k] + CA[k] * CF[k+1]) / (1. + CA[k] + CC[k] - CA[k] * CE[k+1]))
                CFt[k] = (((p0 / grid_p_full[i, j, k])^(Rd_cp) * grid_t[i, j, k] + CA[k] * CFt[k+1]) / (1. + CA[k] + CC[k] - CA[k] .* CE[k+1]))
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
                grid_t[i, j, k] = (CE[k] * grid_t[i, j, k-1] * (p0 / grid_p_full[i, j, k-1])^Rd_cp + CFt[k]) * (grid_p_full[i, j, k] / p0)^Rd_cp
            end
            # --- Backward substitute --- #

        end
    end
    
end
