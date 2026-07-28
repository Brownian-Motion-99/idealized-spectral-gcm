using Test
using JGCM

@testset "Prescribed lower-boundary temperature" begin
    width = 26.0 * pi / 180.0
    expected(latitude) = 29.0 * exp(-latitude^2 / (2.0 * width^2)) + 271.0

    for longitude in (0.0, pi, 1.75pi), latitude in (-0.6, 0.0, 0.4)
        @test Default_Lower_Boundary_Temperature(longitude, latitude) ≈ expected(latitude)
    end

    nλ, nθ, nd = 24, 12, 1
    mesh = Spectral_Spherical_Mesh(3, 4, nλ, nθ, nd, 6.371e6)
    atmo = Atmo_Data(
        "lower-boundary-test",
        nλ,
        nθ,
        nd,
        false,
        false,
        false,
        false,
        mesh.sinθ;
        radius = 6.371e6,
    )

    grid_t = fill(280.0, nλ, nθ, nd)
    grid_shflx = zeros(nλ, nθ, 1)
    wind_speed = fill(10.0, nλ, nθ)
    lowest_level_height = fill(100.0, nλ, nθ)
    density = fill(1.0, nλ, nθ, nd)
    custom_temperature(longitude, latitude) = 300.0 + cos(longitude) + sin(latitude)

    Δt = 10
    C_H = 0.01
    lambda = C_H * wind_speed[1, 1] * Δt / lowest_level_height[1, 1]
    JGCM.Atmos_Param_Module.Sensible_Heating!(
        mesh,
        atmo,
        grid_t,
        grid_shflx,
        wind_speed,
        lowest_level_height,
        density,
        Δt,
        C_H,
        custom_temperature,
    )

    for j = 1:nθ, i = 1:nλ
        expected_t =
            (280.0 + lambda * custom_temperature(mesh.λc[i], mesh.θc[j])) / (1.0 + lambda)
        @test grid_t[i, j, nd] ≈ expected_t
        @test grid_shflx[i, j, 1] ≈
              density[i, j, nd] *
              lowest_level_height[i, j] *
              atmo.cp_air *
              (expected_t - 280.0) / Δt
    end

    grid_ps = fill(100_000.0, nλ, nθ, 1)
    grid_q = fill(0.001, nλ, nθ, nd)
    grid_lhflx = zeros(nλ, nθ, 1)
    C_E = 0.01
    evaporation_lambda = C_E * wind_speed[1, 1] * Δt / lowest_level_height[1, 1]
    JGCM.Atmos_Param_Module.Surface_Evaporation!(
        mesh,
        atmo,
        grid_ps,
        grid_q,
        grid_lhflx,
        wind_speed,
        lowest_level_height,
        density,
        Δt,
        C_E,
        custom_temperature,
    )

    for j = 1:nθ, i = 1:nλ
        surface_temperature = custom_temperature(mesh.λc[i], mesh.θc[j])
        vapor_pressure =
            611.12 * exp(atmo.Lv / atmo.rvgas * (1.0 / 273.15 - 1.0 / surface_temperature))
        surface_q = 0.622 * vapor_pressure / (grid_ps[i, j, 1] - 0.378 * vapor_pressure)
        expected_q = (0.001 + evaporation_lambda * surface_q) / (1.0 + evaporation_lambda)
        @test grid_q[i, j, nd] ≈ expected_q
        @test grid_lhflx[i, j, 1] ≈
              (expected_q - 0.001) *
              density[i, j, nd] *
              lowest_level_height[i, j] *
              atmo.Lv / Δt
    end
end
