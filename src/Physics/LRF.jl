using Base.Threads
using JLD2

abstract type Abstract_LRF_State end

"""Static data for the legacy linear-in-specific-humidity LRF."""
struct LRF_State <: Abstract_LRF_State
    LW_q::Array{Float64,3}   # (nd, nd, nθ)
    ref_q::Array{Float64,3}  # (nλ, nθ, nd)

    function LRF_State(LW_q::AbstractArray{<:Real,3}, ref_q::AbstractArray{<:Real,3})
        nd_out, nd_in, nθ = size(LW_q)
        nλ_ref, nθ_ref, nd_ref = size(ref_q)
        nd_out == nd_in || throw(DimensionMismatch(
            "LRF_LW_q must contain square vertical matrices; got $(size(LW_q))"))
        nθ == nθ_ref || throw(DimensionMismatch(
            "LRF latitude count $nθ does not match ref_q latitude count $nθ_ref"))
        nd_in == nd_ref || throw(DimensionMismatch(
            "LRF vertical size $nd_in does not match ref_q vertical size $nd_ref"))
        nλ_ref > 0 || throw(DimensionMismatch("ref_q must contain at least one longitude"))
        return new(Float64.(LW_q), Float64.(ref_q))
    end
end

"""
    Regularized_LRF_State(B, chi_reference, taper, q0, alpha, nλ)

Static data for `Q = alpha * taper * B * (log(q + q0) - chi_reference)`.
`B` has units K day⁻¹ and `chi_reference` is arranged `(level, latitude)`.
"""
struct Regularized_LRF_State <: Abstract_LRF_State
    B::Matrix{Float64}
    chi_reference::Matrix{Float64}
    taper::Vector{Float64}
    q0::Float64
    alpha::Float64
    nλ::Int
    column_work::Matrix{Float64}

    function Regularized_LRF_State(
        B::AbstractMatrix{<:Real},
        chi_reference::AbstractMatrix{<:Real},
        taper::AbstractVector{<:Real},
        q0::Real,
        alpha::Real,
        nλ::Integer,
    )
        nd_out, nd_in = size(B)
        nd_reference, nθ = size(chi_reference)
        nd_out == nd_in || throw(DimensionMismatch("B_star must be square; got $(size(B))"))
        nd_in == nd_reference || throw(DimensionMismatch(
            "B_star vertical size $nd_in does not match chi_reference size $nd_reference"))
        length(taper) == nθ || throw(DimensionMismatch(
            "taper latitude count $(length(taper)) does not match chi_reference count $nθ"))
        nλ > 0 || throw(DimensionMismatch("LRF grid must contain at least one longitude"))

        B64 = Matrix{Float64}(B)
        reference64 = Matrix{Float64}(chi_reference)
        taper64 = Vector{Float64}(taper)
        all(isfinite, B64) || throw(ArgumentError("B_star must be finite"))
        all(isfinite, reference64) || throw(ArgumentError("chi_reference must be finite"))
        all(x -> isfinite(x) && 0.0 <= x <= 1.0, taper64) ||
            throw(ArgumentError("taper must be finite and lie in [0, 1]"))
        isfinite(q0) && q0 > 0.0 || throw(ArgumentError("q0 must be positive and finite"))
        isfinite(alpha) && alpha >= 0.0 ||
            throw(ArgumentError("alpha must be nonnegative and finite"))

        return new(B64, reference64, taper64, Float64(q0), Float64(alpha), Int(nλ),
                   zeros(Float64, nd_in, nthreads()))
    end
end

"""
    Latitude_LRF_State(B_by_lat, chi_reference, q0, alpha, nλ)

Static data for `Q[:, j] = alpha * B_by_lat[:, :, j] *
(log(q[:, j] + q0) - chi_reference[:, j])`. Each model latitude has its own
vertical response matrix; the same matrix is used at every longitude on that
latitude circle.
"""
struct Latitude_LRF_State <: Abstract_LRF_State
    B_by_lat::Array{Float64,3}
    chi_reference::Matrix{Float64}
    q0::Float64
    alpha::Float64
    nλ::Int
    column_work::Matrix{Float64}

    function Latitude_LRF_State(
        B_by_lat::AbstractArray{<:Real,3},
        chi_reference::AbstractMatrix{<:Real},
        q0::Real,
        alpha::Real,
        nλ::Integer,
    )
        nd_out, nd_in, nθ = size(B_by_lat)
        nd_reference, nθ_reference = size(chi_reference)
        nd_out == nd_in || throw(DimensionMismatch(
            "B_by_lat must contain square vertical matrices; got $(size(B_by_lat))"))
        nd_in == nd_reference || throw(DimensionMismatch(
            "B_by_lat vertical size $nd_in does not match chi_reference size $nd_reference"))
        nθ == nθ_reference || throw(DimensionMismatch(
            "B_by_lat latitude count $nθ does not match chi_reference count $nθ_reference"))
        nλ > 0 || throw(DimensionMismatch("LRF grid must contain at least one longitude"))

        B64 = Array{Float64,3}(B_by_lat)
        reference64 = Matrix{Float64}(chi_reference)
        all(isfinite, B64) || throw(ArgumentError("B_by_lat must be finite"))
        all(isfinite, reference64) || throw(ArgumentError("chi_reference must be finite"))
        isfinite(q0) && q0 > 0.0 || throw(ArgumentError("q0 must be positive and finite"))
        isfinite(alpha) && alpha >= 0.0 ||
            throw(ArgumentError("alpha must be nonnegative and finite"))

        return new(B64, reference64, Float64(q0), Float64(alpha), Int(nλ),
                   zeros(Float64, nd_in, nthreads()))
    end
end

function _require_lrf_key(file, key)
    haskey(file, key) || error("LRF data file is missing variable $key")
    return file[key]
end

function _validate_lrf_latitude(file, latitude, nθ)
    latitude === nothing && return nothing
    stored_latitude = vec(_require_lrf_key(file, "latitude"))
    length(latitude) == nθ || throw(DimensionMismatch(
        "expected $nθ model latitudes, got $(length(latitude))"))
    length(stored_latitude) == nθ || throw(DimensionMismatch(
        "expected stored latitude length $nθ, got $(length(stored_latitude))"))
    all(isapprox.(Float64.(latitude), Float64.(stored_latitude);
                  atol = 1e-10, rtol = 0.0)) ||
        throw(ArgumentError("LRF latitude coordinates do not match the model grid"))
    return nothing
end

"""
    Load_LRF_State(filepath, nλ, nθ, nd; latitude=nothing)

Load a legacy linear LRF, a `regularized_log_tapered_v1` artifact, or a
`regularized_log_latitude_v1` artifact. Supplying model latitudes in degrees
also verifies coordinate values and ordering.
"""
function Load_LRF_State(
    filepath::AbstractString,
    nλ::Int,
    nθ::Int,
    nd::Int;
    latitude = nothing,
)
    isfile(filepath) || error("LRF data file does not exist: $filepath")
    file = JLD2.load(filepath)
    scheme = haskey(file, "scheme") ? String(file["scheme"]) : "linear_q_v1"

    if scheme == "regularized_log_tapered_v1"
        B = _require_lrf_key(file, "B_star")
        chi_reference = _require_lrf_key(file, "chi_reference")
        taper = vec(_require_lrf_key(file, "taper"))
        q0 = _require_lrf_key(file, "q0_kg_kg")
        alpha = _require_lrf_key(file, "alpha")
        size(B) == (nd, nd) || throw(DimensionMismatch(
            "expected B_star size ($nd, $nd), got $(size(B))"))
        size(chi_reference) == (nd, nθ) || throw(DimensionMismatch(
            "expected chi_reference size ($nd, $nθ), got $(size(chi_reference))"))
        length(taper) == nθ || throw(DimensionMismatch(
            "expected taper length $nθ, got $(length(taper))"))

        _validate_lrf_latitude(file, latitude, nθ)
        return Regularized_LRF_State(B, chi_reference, taper, q0, alpha, nλ)
    elseif scheme == "regularized_log_latitude_v1"
        B_by_lat = _require_lrf_key(file, "B_by_lat")
        chi_reference = _require_lrf_key(file, "chi_reference")
        q0 = _require_lrf_key(file, "q0_kg_kg")
        alpha = _require_lrf_key(file, "alpha")
        size(B_by_lat) == (nd, nd, nθ) || throw(DimensionMismatch(
            "expected B_by_lat size ($nd, $nd, $nθ), got $(size(B_by_lat))"))
        size(chi_reference) == (nd, nθ) || throw(DimensionMismatch(
            "expected chi_reference size ($nd, $nθ), got $(size(chi_reference))"))
        _validate_lrf_latitude(file, latitude, nθ)
        return Latitude_LRF_State(B_by_lat, chi_reference, q0, alpha, nλ)
    elseif scheme != "linear_q_v1"
        error("unsupported LRF scheme: $scheme")
    end

    LW_q = _require_lrf_key(file, "LRF_LW_q")
    ref_q = _require_lrf_key(file, "ref_q")
    size(LW_q) == (nd, nd, nθ) || throw(DimensionMismatch(
        "expected LRF_LW_q size ($nd, $nd, $nθ), got $(size(LW_q))"))
    size(ref_q) == (nλ, nθ, nd) || throw(DimensionMismatch(
        "expected ref_q size ($nλ, $nθ, $nd), got $(size(ref_q))"))
    return LRF_State(LW_q, ref_q)
end

function LRF!(
    state::Latitude_LRF_State,
    grid_q::Array{Float64,3},
    grid_lrf_tendency::Array{Float64,3},
    day_to_sec::Int64,
)
    nλ, nθ, nd = size(grid_q)
    expected = (state.nλ, size(state.B_by_lat, 3), size(state.B_by_lat, 1))
    size(grid_q) == expected || throw(DimensionMismatch(
        "grid_q size $(size(grid_q)) does not match LRF grid size $expected"))
    size(grid_lrf_tendency) == size(grid_q) || throw(DimensionMismatch(
        "grid_lrf_tendency size $(size(grid_lrf_tendency)) does not match grid_q size $(size(grid_q))"))
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    invalid_index = findfirst(x -> !isfinite(x) || x < 0.0, grid_q)
    isnothing(invalid_index) || throw(DomainError(
        grid_q[invalid_index],
        "regularized logarithmic LRF requires finite, nonnegative humidity",
    ))

    scale = state.alpha / day_to_sec
    @threads :static for j = 1:nθ
        anomaly = @view state.column_work[:, threadid()]
        for i = 1:nλ
            @inbounds for k = 1:nd
                anomaly[k] = log(grid_q[i, j, k] + state.q0) - state.chi_reference[k, j]
            end
            @inbounds for k_out = 1:nd
                heating_rate = 0.0
                for k_in = 1:nd
                    heating_rate += state.B_by_lat[k_out, k_in, j] * anomaly[k_in]
                end
                grid_lrf_tendency[i, j, k_out] = scale * heating_rate
            end
        end
    end
    return nothing
end

"""Overwrite `grid_lrf_tendency` with the LRF temperature tendency in K s⁻¹."""
function LRF!(
    state::LRF_State,
    grid_q::Array{Float64,3},
    grid_lrf_tendency::Array{Float64,3},
    day_to_sec::Int64,
)
    size(grid_q) == size(state.ref_q) || throw(DimensionMismatch(
        "grid_q size $(size(grid_q)) does not match LRF ref_q size $(size(state.ref_q))"))
    size(grid_lrf_tendency) == size(grid_q) || throw(DimensionMismatch(
        "grid_lrf_tendency size $(size(grid_lrf_tendency)) does not match grid_q size $(size(grid_q))"))
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))

    nλ, nθ, nd = size(grid_q)
    inv_day_to_sec = 1.0 / day_to_sec
    @threads for j = 1:nθ
        for i = 1:nλ
            for k_out = 1:nd
                heating_rate = 0.0
                @inbounds for k_in = 1:nd
                    heating_rate += state.LW_q[k_out, k_in, j] *
                                    (grid_q[i, j, k_in] - state.ref_q[i, j, k_in])
                end
                @inbounds grid_lrf_tendency[i, j, k_out] = heating_rate * inv_day_to_sec
            end
        end
    end
    return nothing
end

function LRF!(
    state::Regularized_LRF_State,
    grid_q::Array{Float64,3},
    grid_lrf_tendency::Array{Float64,3},
    day_to_sec::Int64,
)
    nλ, nθ, nd = size(grid_q)
    expected = (state.nλ, length(state.taper), size(state.B, 1))
    size(grid_q) == expected || throw(DimensionMismatch(
        "grid_q size $(size(grid_q)) does not match LRF grid size $expected"))
    size(grid_lrf_tendency) == size(grid_q) || throw(DimensionMismatch(
        "grid_lrf_tendency size $(size(grid_lrf_tendency)) does not match grid_q size $(size(grid_q))"))
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    invalid_index = findfirst(x -> !isfinite(x) || x < 0.0, grid_q)
    isnothing(invalid_index) || throw(DomainError(
        grid_q[invalid_index],
        "regularized logarithmic LRF requires finite, nonnegative humidity",
    ))

    scale = state.alpha / day_to_sec
    @threads :static for j = 1:nθ
        if iszero(state.taper[j])
            @views fill!(grid_lrf_tendency[:, j, :], 0.0)
            continue
        end
        anomaly = @view state.column_work[:, threadid()]
        latitude_scale = scale * state.taper[j]
        for i = 1:nλ
            @inbounds for k = 1:nd
                anomaly[k] = log(grid_q[i, j, k] + state.q0) - state.chi_reference[k, j]
            end
            @inbounds for k_out = 1:nd
                heating_rate = 0.0
                for k_in = 1:nd
                    heating_rate += state.B[k_out, k_in] * anomaly[k_in]
                end
                grid_lrf_tendency[i, j, k_out] = latitude_scale * heating_rate
            end
        end
    end
    return nothing
end
