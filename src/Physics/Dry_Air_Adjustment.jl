using ...Vert_Coordinate_Module

"""
    Dry_Air_Adjustment!(vert_coord, grid_ps, grid_q_before, grid_q_after,
                        grid_Δp_before)

Adjust moist surface pressure after a fixed-pressure physics update. The
specific humidity at entry to physics (`grid_q_before`) defines the initial dry
air mass, while `grid_q_after` defines the water mass produced by physics.

The model's fixed hybrid coordinate cannot represent the independently
deformed layer thicknesses used by a vertically Lagrangian core. Instead, this
routine changes surface pressure so that dry-air mass is conserved in every
column, then rescales specific humidity to retain the post-physics water mass
in every model layer. All pressure-like mass budgets omit the common factor
`1/g`.
"""
function Dry_Air_Adjustment!(
    vert_coord::Vert_Coordinate,
    grid_ps::Array{Float64,3},
    grid_q_before::Array{Float64,3},
    grid_q_after::Array{Float64,3},
    grid_Δp_before::Array{Float64,3},
)
    nλ, nθ, nd = size(grid_q_after)
    size(grid_ps) == (nλ, nθ, 1) ||
        throw(DimensionMismatch("surface pressure must have size (nλ, nθ, 1)"))
    size(grid_q_before) == (nλ, nθ, nd) ||
        throw(DimensionMismatch("initial humidity does not match adjusted humidity"))
    size(grid_Δp_before) == (nλ, nθ, nd) ||
        throw(DimensionMismatch("initial layer thickness does not match humidity"))
    nd == vert_coord.nd ||
        throw(DimensionMismatch("humidity does not match the vertical coordinate"))

    Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
    surface_pressure_weight = sum(Δbk)
    isfinite(surface_pressure_weight) && surface_pressure_weight > 0.0 ||
        throw(ArgumentError("vertical coordinate must have a positive surface-pressure weight"))

    @inbounds for j in 1:nθ, i in 1:nλ
        ps_before = grid_ps[i, j, 1]
        isfinite(ps_before) && ps_before > 0.0 ||
            throw(DomainError(ps_before, "surface pressure must be finite and positive"))

        dry_before = 0.0
        water_after = 0.0
        water_change = 0.0
        for k in 1:nd
            Δp_before = grid_Δp_before[i, j, k]
            q_before = grid_q_before[i, j, k]
            q_after = grid_q_after[i, j, k]
            isfinite(Δp_before) && Δp_before > 0.0 ||
                throw(DomainError(Δp_before, "layer pressure thickness must be finite and positive"))
            isfinite(q_before) && 0.0 <= q_before < 1.0 ||
                throw(DomainError(q_before, "initial specific humidity must satisfy 0 <= q < 1"))
            isfinite(q_after) && 0.0 <= q_after < 1.0 ||
                throw(DomainError(q_after, "post-physics specific humidity must satisfy 0 <= q < 1"))

            dry_before += (1.0 - q_before) * Δp_before
            water_after += q_after * Δp_before
            water_change += (q_after - q_before) * Δp_before
        end

        ps_after = ps_before + water_change / surface_pressure_weight
        isfinite(ps_after) && ps_after > 0.0 ||
            throw(DomainError(ps_after, "dry-air adjustment produced non-positive surface pressure"))
        grid_ps[i, j, 1] = ps_after

        dry_adjusted = 0.0
        water_adjusted = 0.0
        for k in 1:nd
            Δp_before = grid_Δp_before[i, j, k]
            Δp_after = Δak[k] + Δbk[k] * ps_after
            isfinite(Δp_after) && Δp_after > 0.0 ||
                throw(DomainError(Δp_after, "dry-air adjustment produced a non-positive layer"))

            layer_water = grid_q_after[i, j, k] * Δp_before
            q_adjusted = layer_water / Δp_after
            isfinite(q_adjusted) && 0.0 <= q_adjusted < 1.0 ||
                throw(DomainError(q_adjusted, "adjusted specific humidity must satisfy 0 <= q < 1"))
            grid_q_after[i, j, k] = q_adjusted

            water_adjusted += q_adjusted * Δp_after
            dry_adjusted += (1.0 - q_adjusted) * Δp_after
        end

        budget_scale = max(abs(dry_before), abs(water_after), abs(ps_after), 1.0)
        budget_tolerance = 512.0 * eps(Float64) * budget_scale
        abs(dry_adjusted - dry_before) <= budget_tolerance || error(
            "dry-air adjustment failed dry-mass conservation: residual = " *
            "$(dry_adjusted - dry_before) Pa",
        )
        abs(water_adjusted - water_after) <= budget_tolerance || error(
            "dry-air adjustment failed water conservation: residual = " *
            "$(water_adjusted - water_after) Pa",
        )
    end

    return nothing
end
