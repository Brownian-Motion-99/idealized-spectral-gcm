using Base.Threads
using ...Atmo_Data_Module


function HS_Forcing!(
    atmo_data::Atmo_Data, Δt::Int64, day_to_sec::Int64, sinθ::Array{Float64, 1},
    grid_u::Array{Float64, 3}, grid_v::Array{Float64, 3},
    grid_p_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_t::Array{Float64, 3}, 
    grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3},
    grid_t_eq::Array{Float64, 3}, grid_δt::Array{Float64, 3}, physics_params::Dict{String, Any}
)

    σ_b  = physics_params["σ_b"]::Float64     # Sigma boundary height (unitless)
    k_f  = physics_params["k_f"]::Float64     # Surface friction damping rate (day^{-1})
    k_a  = physics_params["k_a"]::Float64     # Free atmosphere thermal damping rate (day^{-1})
    k_s  = physics_params["k_s"]::Float64     # Surface thermal damping rate (day^{-1})
    ΔT_y = physics_params["ΔT_y"]::Float64    # Equator-to-pole temperature difference (K)
    Δθ_z = physics_params["Δθ_z"]::Float64    # Vertical potential temperature gradient (K)

    cp_air     = atmo_data.cp_air::Float64
    kappa      = atmo_data.kappa::Float64
    nλ, nθ, nd = atmo_data.nλ::Int64, atmo_data.nθ::Int64, atmo_data.nd::Int64

    p_ref = 1.0e5
    t_zero, t_strat  = 315.0, 200.0
    k_a_sec = k_a / day_to_sec
    k_s_sec = k_s / day_to_sec
    k_f_sec = k_f / day_to_sec

    grid_ps = @view grid_p_half[:, :, end]

    inv_one_minus_σb = 1.0 / (1.0 - σ_b)

    @threads for j = 1:nθ
        sinθj = sinθ[j]
        s2    = sinθj^2
        c2    = 1.0 - s2
        c4    = c2^2

        for k = 1:nd
            for i = 1:nλ

                u_val = grid_u[i, j, k]
                v_val = grid_v[i, j, k]
                t_val = grid_t[i, j, k]

                p_s   = grid_ps[i, j]
                p_val = grid_p_full[i, j, k]
                σ     = p_val / p_s

                # --- Rayleigh damping of wind components near the surface --- #
                # Compute Damping Coefficient locally
                # Mathematical logic: max(0, σ - σ_b) handles the condition.
                k_v = 0.0
                if σ > σ_b
                    k_v = k_f_sec * inv_one_minus_σb * (σ - σ_b)
                end

                # Apply Tendencies
                grid_δu[i, j, k] = -k_v * u_val
                grid_δv[i, j, k] = -k_v * v_val
                # --- Rayleigh damping of wind components near the surface --- #

                # --- Frictional heating --- #
                u_term            = (u_val + 0.5 * grid_δu[i, j, k] * Δt) * grid_δu[i, j, k]
                v_term            = (v_val + 0.5 * grid_δv[i, j, k] * Δt) * grid_δv[i, j, k]
                grid_δt[i, j, k] -= (u_term + v_term) / cp_air
                # --- Frictional heating --- #

                # --- Newtonian thermal damping --- #
                # Local Pressure Calculations
                p_norm = p_val / p_ref

                # Compute Equilibrium Temp (t_eq) locally
                term_log           = log(p_norm)
                t_eq_local         = (t_zero - ΔT_y*s2 - Δθ_z*c2*term_log) * (p_norm^kappa)
                t_eq_local         = max(t_strat, t_eq_local)
                grid_t_eq[i, j, k] = t_eq_local

                # Compute Newtonian Damping
                σ   = p_val / p_s
                k_t = k_a_sec
                if σ > σ_b
                    k_t += (k_s_sec - k_a_sec) * inv_one_minus_σb * (σ - σ_b) * c4
                end

                # Update Tendency
                grid_δt[i, j, k] -= k_t * (t_val - t_eq_local)
                # --- Newtonian thermal damping --- #

            end
        end
    end
end
