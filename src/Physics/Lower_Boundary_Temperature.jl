"""
    Default_Lower_Boundary_Temperature(longitude, latitude)

Return the default prescribed lower-boundary temperature in kelvin. This is the
surface temperature used by `Sensible_heat_fluxes!` in the original model with
its zonal perturbation amplitude `A` set to zero.

Both coordinates are in radians. The perturbation is retained explicitly, but
the default result is zonally symmetric because `A = 0`.
"""
function Default_Lower_Boundary_Temperature(longitude::Real, latitude::Real)
    A = 0.0
    σ = 26.0 * pi / 180.0
    σ_damp = 15.0 * pi / 180.0

    T0 = 29.0 * exp(-latitude^2 / (2.0 * σ^2)) + 271.0
    damp = exp(-latitude^2 / (2.0 * σ_damp^2))
    lon_env = 1.0
    P = damp * lon_env * sin(longitude)

    Tsurf = T0 + A * P
    return Tsurf
end
