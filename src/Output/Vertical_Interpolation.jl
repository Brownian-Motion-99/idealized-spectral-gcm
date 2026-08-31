module Vertical_Interpolation_Module

using Base.Threads
using ..Atmo_Data_Module

export Compute_Pressure_Grid!, Interpolate_Field!

const MIN_EXTRAPOLATION_WEIGHT = -0.5
const MAX_EXTRAPOLATION_WEIGHT = 1.5

@inline pressure_log_pressure(p::Float64) = iszero(p) ? 0.0 : p * log(p)

"""
    Compute_Pressure_Grid!(p_full, ak, bk, ps)

Reconstruct full-level pressure using the same hybrid-interface and
Simmons--Burridge definitions used by Isca. `ak` and `bk` are the `nd + 1`
interface coefficients, not arithmetic midpoint coefficients.
"""
function Compute_Pressure_Grid!(
    p_full::AbstractArray{Float64,3},
    ak::AbstractVector{Float64},
    bk::AbstractVector{Float64},
    ps::AbstractArray{Float64,2},
)
    nλ, nθ, nd = size(p_full)
    size(ps) == (nλ, nθ) ||
        throw(DimensionMismatch("surface pressure must match the horizontal pressure grid"))
    length(ak) == nd + 1 && length(bk) == nd + 1 ||
        throw(DimensionMismatch("hybrid interface coefficients must have nd + 1 entries"))
    all(isfinite, ak) && all(isfinite, bk) ||
        throw(ArgumentError("hybrid interface coefficients must be finite"))
    all(p -> isfinite(p) && p > 0.0, ps) ||
        throw(DomainError(ps, "surface pressure must be finite and positive"))

    @threads for j = 1:nθ
        for i = 1:nλ
            @inbounds ps_value = ps[i, j]
            @inbounds for k = 1:nd
                p_top = ak[k] + bk[k] * ps_value
                p_bottom = ak[k+1] + bk[k+1] * ps_value
                Δp = p_bottom - p_top
                if !(p_top >= 0.0 && p_bottom > 0.0 && Δp > 0.0)
                    throw(
                        DomainError(
                            (p_top, p_bottom),
                            "hybrid interface pressures must be non-negative and strictly increasing",
                        ),
                    )
                end

                # This form naturally gives p_bottom/e in a zero-pressure top
                # layer because lim(p*log(p), p -> 0+) = 0.
                log_p_full =
                    (pressure_log_pressure(p_bottom) - pressure_log_pressure(p_top)) / Δp -
                    1.0
                p_full[i, j, k] = exp(log_p_full)
            end
        end
    end
    return nothing
end

"""
    Interpolate_Field!(out, input, p_full, ps, log_targets, var_name, atmo, temperature)

Interpolate a model-level field linearly in log pressure following Isca's
pressure-level postprocessor. Targets below the lowest model full level are
masked with `NaN`; limited extrapolation is retained above the top model level.
Geopotential height uses Isca's hydrostatic interpolation when temperature is
available.

The `ps` argument is retained for API compatibility and to document that
`p_full` was reconstructed from the surface pressure for the same output
record. Temporal averaging is deliberately performed by `Output_Manager`
before this routine is called.
"""
function Interpolate_Field!(
    out_3d::AbstractArray{Float64,3},
    in_3d::AbstractArray{Float64,3},
    p_3d::AbstractArray{Float64,3},
    ps_2d::AbstractArray{Float64,2},
    log_targets::AbstractVector{Float64},
    var_name::Symbol,
    atmo_data::Atmo_Data,
    t_3d::Union{AbstractArray{Float64,3},Nothing} = nothing,
)
    nλ, nθ, n_plev = size(out_3d)
    nd = size(in_3d, 3)
    size(in_3d) == size(p_3d) ||
        throw(DimensionMismatch("field and full-level pressure grids must match"))
    size(in_3d, 1) == nλ && size(in_3d, 2) == nθ ||
        throw(DimensionMismatch("input and output horizontal grids must match"))
    size(ps_2d) == (nλ, nθ) ||
        throw(DimensionMismatch("surface pressure must match the horizontal grid"))
    length(log_targets) == n_plev ||
        throw(DimensionMismatch("target pressure count must match output levels"))
    if var_name == :z
        isnothing(t_3d) &&
            throw(ArgumentError("temperature is required to interpolate height"))
        size(t_3d) == size(in_3d) ||
            throw(DimensionMismatch("temperature and height grids must match"))
    end

    rd_over_grav = atmo_data.rdgas / atmo_data.grav

    @threads for j = 1:nθ
        for i = 1:nλ
            @inbounds begin
                if nd == 1
                    for n = 1:n_plev
                        out_3d[i, j, n] =
                            log_targets[n] <= log(p_3d[i, j, 1]) ? in_3d[i, j, 1] : NaN
                    end
                    continue
                end

                for n = 1:n_plev
                    log_target = log_targets[n]

                    # As in Isca's default plevel.sh path, pressure levels below
                    # the lowest model full level are treated as missing.
                    if log_target > log(p_3d[i, j, nd])
                        out_3d[i, j, n] = NaN
                        continue
                    end

                    k_bottom = 2
                    for k = 2:nd
                        if log_target <= log(p_3d[i, j, k])
                            k_bottom = k
                            break
                        end
                    end

                    log_p_bottom = log(p_3d[i, j, k_bottom])
                    log_p_top = log(p_3d[i, j, k_bottom-1])
                    weight = (log_target - log_p_bottom) / (log_p_top - log_p_bottom)
                    weight =
                        clamp(weight, MIN_EXTRAPOLATION_WEIGHT, MAX_EXTRAPOLATION_WEIGHT)

                    value_bottom = in_3d[i, j, k_bottom]
                    value_top = in_3d[i, j, k_bottom-1]
                    interpolated = value_bottom + weight * (value_top - value_bottom)

                    if var_name == :z
                        temperature_bottom = t_3d[i, j, k_bottom]
                        temperature_top = t_3d[i, j, k_bottom-1]
                        temperature_target =
                            temperature_bottom +
                            weight * (temperature_top - temperature_bottom)
                        out_3d[i, j, n] =
                            value_bottom +
                            (log_p_bottom - log_target) *
                            (temperature_target + temperature_bottom) *
                            (0.5 * rd_over_grav)
                    else
                        out_3d[i, j, n] = interpolated
                    end
                end
            end
        end
    end
    return nothing
end

end
