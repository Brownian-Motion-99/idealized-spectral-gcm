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
    grid_p_half = zeros(nλ, nθ, nd + 1)
    grid_p_half[:, :, 1] .= 80_000.0
    grid_p_half[:, :, 2] .= 100_000.0
    custom_temperature(longitude, latitude) = 300.0 + cos(longitude) + sin(latitude)

    Δt = 10
    C_H = 0.01
    lambda = C_H * wind_speed[1, 1] * Δt / lowest_level_height[1, 1]
    JGCM.Atmos_Param_Module.Sensible_Heating!(
        mesh,
        atmo,
        grid_p_half,
        grid_t,
        grid_shflx,
        wind_speed,
        lowest_level_height,
        Δt,
        C_H,
        custom_temperature,
    )

    expected_t = similar(grid_t)
    for j = 1:nθ, i = 1:nλ
        expected_t[i, j, nd] =
            (280.0 + lambda * custom_temperature(mesh.λc[i], mesh.θc[j])) / (1.0 + lambda)
    end
    expected_shflx =
        @. ((grid_p_half[:, :, nd+1] - grid_p_half[:, :, nd]) / atmo.grav) * atmo.cp_air *
           (expected_t[:, :, nd] - 280.0) / Δt
    @test grid_t ≈ expected_t
    @test grid_shflx[:, :, 1] ≈ expected_shflx

    grid_ps = fill(100_000.0, nλ, nθ, 1)
    grid_q = fill(0.001, nλ, nθ, nd)
    grid_lhflx = zeros(nλ, nθ, 1)
    C_E = 0.01
    evaporation_lambda = C_E * wind_speed[1, 1] * Δt / lowest_level_height[1, 1]
    JGCM.Atmos_Param_Module.Surface_Evaporation!(
        mesh,
        atmo,
        grid_ps,
        grid_p_half,
        grid_q,
        grid_lhflx,
        wind_speed,
        lowest_level_height,
        Δt,
        C_E,
        custom_temperature,
    )

    expected_q = similar(grid_q)
    for j = 1:nθ, i = 1:nλ
        surface_temperature = custom_temperature(mesh.λc[i], mesh.θc[j])
        surface_q = Saturation_Specific_Humidity(
            surface_temperature,
            grid_ps[i, j, 1],
            atmo.rdgas / atmo.rvgas,
        )
        expected_q[i, j, nd] =
            (0.001 + evaporation_lambda * surface_q) / (1.0 + evaporation_lambda)
    end
    expected_lhflx =
        @. (expected_q[:, :, nd] - 0.001) *
           ((grid_p_half[:, :, nd+1] - grid_p_half[:, :, nd]) / atmo.grav) * atmo.Lv / Δt
    @test grid_q ≈ expected_q
    @test grid_lhflx[:, :, 1] ≈ expected_lhflx
end
