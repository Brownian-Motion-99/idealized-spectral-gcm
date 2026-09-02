module Vertical_Interpolation_Module

using Base.Threads
using ..Atmo_Data_Module

export Pressure_Interpolation_Cache, Prepare_Interpolation!, Apply_Interpolation!
export Compute_Pressure_Grid!, Interpolate_Field!

const MIN_EXTRAPOLATION_WEIGHT = -0.5
const MAX_EXTRAPOLATION_WEIGHT = 1.5

@inline pressure_log_pressure(p::Float64) = iszero(p) ? 0.0 : p * log(p)

"""
    Pressure_Interpolation_Cache(nlon, nlat, nlevels, ntargets)

Reusable interpolation geometry for one pressure field. As in Isca's
`pres_interp_type`, the expensive logarithms, vertical searches, and weights
are prepared once and then shared by every field in an output record.
"""
struct Pressure_Interpolation_Cache
    log_p_full::Array{Float64,3}
    log_targets::Vector{Float64}
    k_bottom::Array{Int,3}
    weight::Array{Float64,3}
end

function Pressure_Interpolation_Cache(nλ::Integer, nθ::Integer, nd::Integer, n_plev::Integer)
    nλ >= 0 && nθ >= 0 && nd >= 1 && n_plev >= 0 ||
        throw(ArgumentError("interpolation-cache dimensions must be non-negative and nd positive"))
    return Pressure_Interpolation_Cache(
        zeros(Float64, nλ, nθ, nd),
        zeros(Float64, n_plev),
        zeros(Int, nλ, nθ, n_plev),
        zeros(Float64, nλ, nθ, n_plev),
    )
end

function validate_hybrid_grid(ak, bk, ps, nλ, nθ, nd)
    size(ps) == (nλ, nθ) ||
        throw(DimensionMismatch("surface pressure must match the horizontal pressure grid"))
    length(ak) == nd + 1 && length(bk) == nd + 1 ||
        throw(DimensionMismatch("hybrid interface coefficients must have nd + 1 entries"))
    all(isfinite, ak) && all(isfinite, bk) ||
        throw(ArgumentError("hybrid interface coefficients must be finite"))
    all(p -> isfinite(p) && p > 0.0, ps) ||
        throw(DomainError(ps, "surface pressure must be finite and positive"))
    return nothing
end

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
    validate_hybrid_grid(ak, bk, ps, nλ, nθ, nd)

    @threads for j = 1:nθ
        for i = 1:nλ
            @inbounds ps_value = ps[i, j]
            p_top = ak[1] + bk[1] * ps_value
            p_top >= 0.0 || throw(
                DomainError(p_top, "hybrid interface pressures must be non-negative"),
            )
            p_top_log = pressure_log_pressure(p_top)
            @inbounds for k = 1:nd
                p_bottom = ak[k+1] + bk[k+1] * ps_value
                Δp = p_bottom - p_top
                p_bottom > 0.0 && Δp > 0.0 || throw(
                    DomainError(
                        (p_top, p_bottom),
                        "hybrid interface pressures must be non-negative and strictly increasing",
                    ),
                )
                p_bottom_log = pressure_log_pressure(p_bottom)
                p_full[i, j, k] = exp((p_bottom_log - p_top_log) / Δp - 1.0)
                p_top = p_bottom
                p_top_log = p_bottom_log
            end
        end
    end
    return nothing
end

function validate_cache(cache::Pressure_Interpolation_Cache, nλ, nθ, nd, n_plev)
    size(cache.log_p_full) == (nλ, nθ, nd) ||
        throw(DimensionMismatch("cached full-level pressure grid has the wrong dimensions"))
    length(cache.log_targets) == n_plev ||
        throw(DimensionMismatch("cached pressure targets have the wrong dimensions"))
    size(cache.k_bottom) == (nλ, nθ, n_plev) ||
        throw(DimensionMismatch("cached interpolation indices have the wrong dimensions"))
    size(cache.weight) == (nλ, nθ, n_plev) ||
        throw(DimensionMismatch("cached interpolation weights have the wrong dimensions"))
    return nothing
end

function prepare_indices_and_weights!(cache::Pressure_Interpolation_Cache, log_targets)
    nλ, nθ, nd = size(cache.log_p_full)
    n_plev = length(log_targets)
    validate_cache(cache, nλ, nθ, nd, n_plev)
    all(isfinite, log_targets) || throw(ArgumentError("log pressure targets must be finite"))
    copyto!(cache.log_targets, log_targets)

    @threads for j = 1:nθ
        @inbounds for n = 1:n_plev
            log_target = cache.log_targets[n]
            for i = 1:nλ
                if log_target > cache.log_p_full[i, j, nd]
                    # Zero marks pressure levels below the lowest model full
                    # level, which are masked in Isca's default path.
                    cache.k_bottom[i, j, n] = 0
                    cache.weight[i, j, n] = 0.0
                    continue
                end

                if nd == 1
                    cache.k_bottom[i, j, n] = 1
                    cache.weight[i, j, n] = 0.0
                    continue
                end

                k_bottom = 2
                for k = 2:nd
                    if log_target <= cache.log_p_full[i, j, k]
                        k_bottom = k
                        break
                    end
                end
                log_p_bottom = cache.log_p_full[i, j, k_bottom]
                log_p_top = cache.log_p_full[i, j, k_bottom-1]
                cache.k_bottom[i, j, n] = k_bottom
                cache.weight[i, j, n] = clamp(
                    (log_target - log_p_bottom) / (log_p_top - log_p_bottom),
                    MIN_EXTRAPOLATION_WEIGHT,
                    MAX_EXTRAPOLATION_WEIGHT,
                )
            end
        end
    end
    return nothing
end

"""
    Prepare_Interpolation!(cache, ak, bk, ps, log_targets)

Prepare the pressure geometry for an output record directly from the hybrid
coordinate. This avoids materializing pressure and then taking its logarithm.
Call it once per surface-pressure record, then use `Apply_Interpolation!` for
every three-dimensional variable.
"""
function Prepare_Interpolation!(
    cache::Pressure_Interpolation_Cache,
    ak::AbstractVector{Float64},
    bk::AbstractVector{Float64},
    ps::AbstractArray{Float64,2},
    log_targets::AbstractVector{Float64},
)
    nλ, nθ, nd = size(cache.log_p_full)
    validate_hybrid_grid(ak, bk, ps, nλ, nθ, nd)
    validate_cache(cache, nλ, nθ, nd, length(log_targets))

    @threads for j = 1:nθ
        for i = 1:nλ
            @inbounds ps_value = ps[i, j]
            p_top = ak[1] + bk[1] * ps_value
            p_top >= 0.0 || throw(
                DomainError(p_top, "hybrid interface pressures must be non-negative"),
            )
            p_top_log = pressure_log_pressure(p_top)
            @inbounds for k = 1:nd
                p_bottom = ak[k+1] + bk[k+1] * ps_value
                Δp = p_bottom - p_top
                p_bottom > 0.0 && Δp > 0.0 || throw(
                    DomainError(
                        (p_top, p_bottom),
                        "hybrid interface pressures must be non-negative and strictly increasing",
                    ),
                )
                p_bottom_log = pressure_log_pressure(p_bottom)
                cache.log_p_full[i, j, k] = (p_bottom_log - p_top_log) / Δp - 1.0
                p_top = p_bottom
                p_top_log = p_bottom_log
            end
        end
    end
    prepare_indices_and_weights!(cache, log_targets)
    return nothing
end

"""
    Prepare_Interpolation!(cache, p_full, log_targets)

Prepare interpolation geometry from an existing full-level pressure array.
This overload supports callers which already have pressure in memory.
"""
function Prepare_Interpolation!(
    cache::Pressure_Interpolation_Cache,
    p_full::AbstractArray{Float64,3},
    log_targets::AbstractVector{Float64},
)
    nλ, nθ, nd = size(p_full)
    validate_cache(cache, nλ, nθ, nd, length(log_targets))
    all(p -> isfinite(p) && p > 0.0, p_full) ||
        throw(DomainError(p_full, "full-level pressures must be finite and positive"))

    @threads for j = 1:nθ
        @inbounds for k = 1:nd, i = 1:nλ
            cache.log_p_full[i, j, k] = log(p_full[i, j, k])
        end
    end
    prepare_indices_and_weights!(cache, log_targets)
    return nothing
end

"""
    Apply_Interpolation!(out, input, cache, var_name, atmo, temperature=nothing)

Apply precomputed pressure interpolation geometry to one field. Geopotential
height uses Isca's hydrostatic interpolation when temperature is supplied.
"""
function Apply_Interpolation!(
    out_3d::AbstractArray{Float64,3},
    in_3d::AbstractArray{Float64,3},
    cache::Pressure_Interpolation_Cache,
    var_name::Symbol,
    atmo_data::Atmo_Data,
    t_3d::Union{AbstractArray{Float64,3},Nothing} = nothing,
)
    nλ, nθ, n_plev = size(out_3d)
    nd = size(in_3d, 3)
    size(in_3d, 1) == nλ && size(in_3d, 2) == nθ ||
        throw(DimensionMismatch("input and output horizontal grids must match"))
    validate_cache(cache, nλ, nθ, nd, n_plev)
    if var_name == :z
        isnothing(t_3d) && throw(ArgumentError("temperature is required to interpolate height"))
        size(t_3d) == size(in_3d) ||
            throw(DimensionMismatch("temperature and height grids must match"))
    end

    half_rd_over_grav = 0.5 * atmo_data.rdgas / atmo_data.grav

    # Longitude is innermost to follow Julia's column-major memory layout for
    # the output and cache arrays.
    @threads for j = 1:nθ
        @inbounds for n = 1:n_plev
            log_target = cache.log_targets[n]
            for i = 1:nλ
                k_bottom = cache.k_bottom[i, j, n]
                if iszero(k_bottom)
                    out_3d[i, j, n] = NaN
                    continue
                elseif nd == 1
                    out_3d[i, j, n] = in_3d[i, j, 1]
                    continue
                end

                weight = cache.weight[i, j, n]
                value_bottom = in_3d[i, j, k_bottom]
                value_top = in_3d[i, j, k_bottom-1]
                if var_name == :z
                    temperature_bottom = t_3d[i, j, k_bottom]
                    temperature_top = t_3d[i, j, k_bottom-1]
                    temperature_target =
                        temperature_bottom + weight * (temperature_top - temperature_bottom)
                    out_3d[i, j, n] =
                        value_bottom +
                        (cache.log_p_full[i, j, k_bottom] - log_target) *
                        (temperature_target + temperature_bottom) * half_rd_over_grav
                else
                    out_3d[i, j, n] = value_bottom + weight * (value_top - value_bottom)
                end
            end
        end
    end
    return nothing
end

"""
    Interpolate_Field!(out, input, p_full, ps, log_targets, var_name, atmo, temperature)

Compatibility wrapper for one-off callers. Hot paths should prepare a reusable
`Pressure_Interpolation_Cache` once and call `Apply_Interpolation!` for each
field. Temporal averaging remains the responsibility of `Output_Manager`.
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
    size(p_3d) == size(in_3d) ||
        throw(DimensionMismatch("field and full-level pressure grids must match"))
    size(ps_2d) == (nλ, nθ) ||
        throw(DimensionMismatch("surface pressure must match the horizontal grid"))
    length(log_targets) == n_plev ||
        throw(DimensionMismatch("target pressure count must match output levels"))

    cache = Pressure_Interpolation_Cache(nλ, nθ, size(in_3d, 3), n_plev)
    Prepare_Interpolation!(cache, p_3d, log_targets)
    Apply_Interpolation!(out_3d, in_3d, cache, var_name, atmo_data, t_3d)
    return nothing
end

end
