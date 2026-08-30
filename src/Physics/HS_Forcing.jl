using Base.Threads
using ...Atmo_Data_Module



"""
    Rayleigh_Friction!(atmo_data, effective_dt, day_to_sec,
                       p_half, p_full, u, v, temperature, physics_params)

Apply the low-level Held--Suarez Rayleigh damping as a finite forward-Euler
step on a physics working state. The kinetic energy removed by damping is
returned locally as frictional heating, so the combined wind/temperature
update conserves kinetic plus sensible energy to roundoff.
"""
function Rayleigh_Friction!(
    atmo_data::Atmo_Data,
    effective_dt::Real,
    day_to_sec::Integer,
    grid_p_half::Array{Float64,3},
    grid_p_full::Array{Float64,3},
    grid_u::Array{Float64,3},
    grid_v::Array{Float64,3},
    grid_t::Array{Float64,3},
    physics_params::Dict{String,Any},
)
    size(grid_u) == size(grid_v) == size(grid_t) == size(grid_p_full) ||
        throw(DimensionMismatch("Rayleigh-friction full-level fields must match"))
    nλ, nθ, nd = size(grid_u)
    size(grid_p_half) == (nλ, nθ, nd + 1) ||
        throw(DimensionMismatch("Rayleigh-friction p_half has incorrect size"))

    dt = Float64(effective_dt)
    isfinite(dt) && dt > 0 || throw(ArgumentError("effective_dt must be positive and finite"))
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    sigma_b = Float64(physics_params["σ_b"])
    k_f = Float64(physics_params["k_f"]) / Float64(day_to_sec)
    0.0 <= sigma_b < 1.0 || throw(ArgumentError("σ_b must lie in [0, 1)"))
    isfinite(k_f) && k_f >= 0 || throw(ArgumentError("k_f must be finite and non-negative"))
    cp = atmo_data.cp_air
    inv_one_minus_sigma_b = 1.0 / (1.0 - sigma_b)
    grid_ps = @view grid_p_half[:, :, end]

    @threads for j = 1:nθ
        for k = 1:nd
            for i = 1:nλ
                sigma = grid_p_full[i, j, k] / grid_ps[i, j]
                damping_rate =
                    sigma > sigma_b ? k_f * (sigma - sigma_b) * inv_one_minus_sigma_b : 0.0
                damping_fraction = dt * damping_rate
                damping_fraction <= 1.0 || throw(
                    ArgumentError(
                        "explicit Rayleigh damping requires k_v * effective_dt <= 1; " *
                        "got $damping_fraction at (i=$i, j=$j, k=$k)",
                    ),
                )

                u_old = grid_u[i, j, k]
                v_old = grid_v[i, j, k]
                u_new = (1.0 - damping_fraction) * u_old
                v_new = (1.0 - damping_fraction) * v_old
                grid_u[i, j, k] = u_new
                grid_v[i, j, k] = v_new
                grid_t[i, j, k] +=
                    0.5 * (u_old^2 + v_old^2 - u_new^2 - v_new^2) / cp
            end
        end
    end
    return nothing
end

"""
    Newtonian_Relaxation!(atmo_data, effective_dt, day_to_sec, sinθ,
                           p_half, p_full, temperature, t_eq, physics_params)

Apply the modified Held--Suarez Newtonian relaxation to the post-friction
physics working temperature. `t_eq` is overwritten as a diagnostic.
"""
function Newtonian_Relaxation!(
    atmo_data::Atmo_Data,
    effective_dt::Real,
    day_to_sec::Integer,
    sinθ::AbstractVector{<:Real},
    grid_p_half::Array{Float64,3},
    grid_p_full::Array{Float64,3},
    grid_t::Array{Float64,3},
    grid_t_eq::Array{Float64,3},
    physics_params::Dict{String,Any},
)
    size(grid_t) == size(grid_p_full) == size(grid_t_eq) ||
        throw(DimensionMismatch("Newtonian-relaxation full-level fields must match"))
    nλ, nθ, nd = size(grid_t)
    length(sinθ) == nθ || throw(DimensionMismatch("sinθ has incorrect length"))
    size(grid_p_half) == (nλ, nθ, nd + 1) ||
        throw(DimensionMismatch("Newtonian-relaxation p_half has incorrect size"))

    dt = Float64(effective_dt)
    isfinite(dt) && dt > 0 || throw(ArgumentError("effective_dt must be positive and finite"))
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    sigma_b = Float64(physics_params["σ_b"])
    k_a = Float64(physics_params["k_a"]) / Float64(day_to_sec)
    k_s = Float64(physics_params["k_s"]) / Float64(day_to_sec)
    delta_t_y = Float64(physics_params["ΔT_y"])
    delta_theta_z = Float64(physics_params["Δθ_z"])
    t_equator = Float64(get(physics_params, "T_equator", 315.0))
    t_stratosphere = Float64(get(physics_params, "T_stratosphere", 200.0))
    0.0 <= sigma_b < 1.0 || throw(ArgumentError("σ_b must lie in [0, 1)"))
    all(x -> isfinite(x) && x >= 0, (k_a, k_s)) ||
        throw(ArgumentError("HS thermal damping rates must be finite and non-negative"))
    isfinite(t_equator) && t_equator > 0 ||
        throw(ArgumentError("T_equator must be positive and finite"))
    isfinite(t_stratosphere) && t_stratosphere > 0 ||
        throw(ArgumentError("T_stratosphere must be positive and finite"))

    kappa = atmo_data.kappa
    p_ref = 1.0e5
    inv_one_minus_sigma_b = 1.0 / (1.0 - sigma_b)
    grid_ps = @view grid_p_half[:, :, end]

    @threads for j = 1:nθ
        sin_lat = Float64(sinθ[j])
        sin2 = sin_lat^2
        cos2 = 1.0 - sin2
        cos4 = cos2^2
        for k = 1:nd
            for i = 1:nλ
                pressure = grid_p_full[i, j, k]
                p_norm = pressure / p_ref
                p_norm > 0 || throw(ArgumentError("full-level pressure must be positive"))
                t_eq =
                    (t_equator - delta_t_y * sin2 - delta_theta_z * cos2 * log(p_norm)) *
                    p_norm^kappa
                t_eq = max(t_stratosphere, t_eq)
                grid_t_eq[i, j, k] = t_eq

                sigma = pressure / grid_ps[i, j]
                damping_rate = k_a
                if sigma > sigma_b
                    damping_rate +=
                        (k_s - k_a) * (sigma - sigma_b) * inv_one_minus_sigma_b * cos4
                end
                relaxation_fraction = dt * damping_rate
                0.0 <= relaxation_fraction <= 1.0 || throw(
                    ArgumentError(
                        "explicit Newtonian relaxation requires 0 <= k_t * effective_dt <= 1; " *
                        "got $relaxation_fraction at (i=$i, j=$j, k=$k)",
                    ),
                )
                grid_t[i, j, k] += relaxation_fraction * (t_eq - grid_t[i, j, k])
            end
        end
    end
    return nothing
end
