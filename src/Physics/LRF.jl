using Base.Threads
using JLD2

"""
    LRF_State

Static data used by the water vapor linear-response forcing. `LW_q` contains one
`nd × nd` response matrix per latitude, and `ref_q` is the reference specific
water vapor specific humidity on the model grid.
"""
struct LRF_State
    LW_q::Array{Float64,3}   # (nd, nd, nθ)
    ref_q::Array{Float64,3}  # (nλ, nθ, nd)

    function LRF_State(LW_q::AbstractArray{<:Real,3}, ref_q::AbstractArray{<:Real,3})
        nd_out, nd_in, nθ = size(LW_q)
        nλ_ref, nθ_ref, nd_ref = size(ref_q)

        nd_out == nd_in || throw(
            DimensionMismatch("LRF_LW_q must contain square vertical matrices; got $(size(LW_q))"),
        )
        nθ == nθ_ref || throw(
            DimensionMismatch(
                "LRF latitude count $nθ does not match ref_q latitude count $nθ_ref",
            ),
        )
        nd_in == nd_ref || throw(
            DimensionMismatch(
                "LRF vertical size $nd_in does not match ref_q vertical size $nd_ref",
            ),
        )
        nλ_ref > 0 || throw(DimensionMismatch("ref_q must contain at least one longitude"))

        return new(Float64.(LW_q), Float64.(ref_q))
    end
end

"""
    Load_LRF_State(filepath, nλ, nθ, nd)

Load and validate the time-independent LRF data once, before time integration.
The legacy data layout is `(nd, nd, nθ)` for `LRF_LW_q` and `(nλ, nθ, nd)`
for `ref_q`.
"""
function Load_LRF_State(filepath::AbstractString, nλ::Int, nθ::Int, nd::Int)
    isfile(filepath) || error("LRF data file does not exist: $filepath")

    file = JLD2.load(filepath)
    haskey(file, "LRF_LW_q") || error("LRF data file is missing variable LRF_LW_q")
    haskey(file, "ref_q") || error("LRF data file is missing variable ref_q")

    LW_q = file["LRF_LW_q"]
    ref_q = file["ref_q"]

    size(LW_q) == (nd, nd, nθ) || throw(
        DimensionMismatch(
            "expected LRF_LW_q size ($nd, $nd, $nθ), got $(size(LW_q))",
        ),
    )
    size(ref_q) == (nλ, nθ, nd) || throw(
        DimensionMismatch(
            "expected ref_q size ($nλ, $nθ, $nd), got $(size(ref_q))",
        ),
    )

    return LRF_State(LW_q, ref_q)
end

"""
    LRF!(state, grid_q, grid_lrf_tendency, day_to_sec)

Calculate the columnwise temperature tendency caused by water vapor anomalies.
`grid_lrf_tendency` is overwritten in-place and returned in K s⁻¹.
"""
function LRF!(
    state::LRF_State,
    grid_q::Array{Float64,3},
    grid_lrf_tendency::Array{Float64,3},
    day_to_sec::Int64,
)
    size(grid_q) == size(state.ref_q) || throw(
        DimensionMismatch(
            "grid_q size $(size(grid_q)) does not match LRF ref_q size $(size(state.ref_q))",
        ),
    )
    size(grid_lrf_tendency) == size(grid_q) || throw(
        DimensionMismatch(
            "grid_lrf_tendency size $(size(grid_lrf_tendency)) does not match grid_q size $(size(grid_q))",
        ),
    )
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))

    nλ, nθ, nd = size(grid_q)
    inv_day_to_sec = 1.0 / day_to_sec

    @threads for j = 1:nθ
        for i = 1:nλ
            for k_out = 1:nd
                heating_rate = 0.0
                @inbounds for k_in = 1:nd
                    heating_rate +=
                        state.LW_q[k_out, k_in, j] *
                        (grid_q[i, j, k_in] - state.ref_q[i, j, k_in])
                end
                @inbounds grid_lrf_tendency[i, j, k_out] =
                    heating_rate * inv_day_to_sec
            end
        end
    end

    return nothing
end
