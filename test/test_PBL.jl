@testset "PBL surface-layer geometry and density" begin
    nλ, nθ, nd = 2, 2, 2
    atmo = Atmo_Data(
        "PBL test", nλ, nθ, nd,
        false, false, false, true,
        zeros(nθ);
        radius = 6.371e6,
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
    workspace = JGCM.Atmos_Param_Module.PBL_Workspace(nλ, nθ, nd)

    wind, za, rho = JGCM.Atmos_Param_Module.Calculate_V_c_za_rho!(
        workspace, atmo, p_half, p_full, ps, u, v, temperature, humidity,
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

@testset "PBL mixing uses post-state interface density" begin
    nλ, nθ, nd = 1, 1, 2
    atmo = Atmo_Data(
        "PBL interface-density test", nλ, nθ, nd,
        false, false, false, true,
        zeros(nθ);
        radius = 6.371e6,
    )
    p_half = reshape([0.0, 50_000.0, 100_000.0], nλ, nθ, nd + 1)
    p_full = reshape([25_000.0, 75_000.0], nλ, nθ, nd)
    theta = 300.0
    temperature = reshape(
        [theta * (p_full[k] / 100_000.0)^atmo.kappa for k = 1:nd],
        nλ,
        nθ,
        nd,
    )
    humidity = reshape([0.001, 0.010], nλ, nθ, nd)
    humidity_initial = copy(humidity)
    K_E = zeros(nλ, nθ, nd + 1)
    surface_wind = fill(10.0, nλ, nθ)
    surface_height = fill(100.0, nλ, nθ)
    effective_dt = 600
    coefficient = 0.0044
    params = Dict{String,Any}(
        "PBL_Top_Mode" => :PressureLevel,
        "PBL_Top_Value" => 0.0,
    )

    tv_upper = temperature[1] * (1.0 + 0.608 * humidity[1])
    tv_lower = temperature[2] * (1.0 + 0.608 * humidity[2])
    rho_interface = p_half[2] / (atmo.rdgas * 0.5 * (tv_upper + tv_lower))
    diffusivity = coefficient * surface_wind[1] * surface_height[1]
    ca =
        effective_dt * atmo.grav^2 * diffusivity * rho_interface^2 /
        ((p_half[2] - p_half[1]) * (p_full[2] - p_full[1]))
    cc =
        effective_dt * atmo.grav^2 * diffusivity * rho_interface^2 /
        ((p_half[3] - p_half[2]) * (p_full[2] - p_full[1]))
    denominator = 1.0 + ca + cc
    expected_q_upper =
        ((1.0 + cc) * humidity_initial[1] + ca * humidity_initial[2]) / denominator
    expected_q_lower =
        (cc * humidity_initial[1] + (1.0 + ca) * humidity_initial[2]) / denominator

    JGCM.Atmos_Param_Module.Implicit_PBL_Mixing!(
        atmo,
        p_full,
        p_half,
        temperature,
        humidity,
        K_E,
        surface_wind,
        surface_height,
        params,
        effective_dt,
        coefficient,
    )

    @test humidity[1] ≈ expected_q_upper
    @test humidity[2] ≈ expected_q_lower
    @test humidity[1] > humidity_initial[1]
    @test humidity[2] < humidity_initial[2]
    @test temperature ≈
          reshape(
        [theta * (p_full[k] / 100_000.0)^atmo.kappa for k = 1:nd],
        nλ,
        nθ,
        nd,
    )
end
