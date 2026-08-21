@testset "PBL surface-layer geometry and density" begin
    nλ, nθ, nd = 2, 2, 2
    atmo = Atmo_Data(
        "PBL test", nλ, nθ, nd,
        false, false, false, true,
        zeros(nθ);
        radius=6.371e6,
    )

    # Ordering is top to bottom. The upper interface of the lowest layer is
    # 90 kPa and the surface is 100 kPa.
    p_half = zeros(nλ, nθ, nd + 1)
    p_half[:, :, 1] .= 40_000.0
    p_half[:, :, 2] .= 90_000.0
    p_half[:, :, 3] .= 100_000.0

    p_full = zeros(nλ, nθ, nd)
    p_full[:, :, 1] .= 60_000.0
    p_full[:, :, 2] .= 95_000.0

    ps = fill(100_000.0, nλ, nθ, 1)
    u = fill(3.0, nλ, nθ, nd)
    v = fill(4.0, nλ, nθ, nd)
    temperature = zeros(nλ, nθ, nd)
    temperature[:, :, 1] .= 250.0
    temperature[:, :, 2] .= 280.0
    humidity = zeros(nλ, nθ, nd)
    humidity[:, :, 1] .= 0.001
    humidity[:, :, 2] .= 0.010

    wind, za, rho = JGCM.Atmos_Param_Module.Calculate_V_c_za_rho(
        atmo, p_half, p_full, ps, u, v, temperature, humidity,
    )

    expected_za = (atmo.rdgas / atmo.grav) * 280.0 * (1.0 + 0.608 * 0.010) *
                  log(100_000.0 / 90_000.0) / 2.0
    expected_rho_upper = 60_000.0 / (atmo.rdgas * 250.0 * (1.0 + 0.608 * 0.001))
    expected_rho_lower = 95_000.0 / (atmo.rdgas * 280.0 * (1.0 + 0.608 * 0.010))

    @test all(wind .== 5.0)
    @test all(za .≈ expected_za)
    @test all(rho[:, :, 1] .≈ expected_rho_upper)
    @test all(rho[:, :, 2] .≈ expected_rho_lower)
end
