using JGCM

function lscale_test_atmosphere(nλ, nθ, nd)
    return Atmo_Data(
        "lscale_test",
        nλ,
        nθ,
        nd,
        false,
        false,
        false,
        false,
        zeros(nθ);
        radius = 6.371e6,
    )
end

function lscale_saturation_specific_humidity(atmo, temperature, pressure)
    epsilon = atmo.rdgas / atmo.rvgas
    mixing_ratio = epsilon * Betts_Miller_Saturation_Vapor_Pressure(temperature) / pressure
    return mixing_ratio / (1.0 + mixing_ratio)
end

function run_lscale_column(effective_dt, heating_scale)
    nλ, nθ, nd = 1, 1, 2
    atmo = lscale_test_atmosphere(nλ, nθ, nd)
    p_half = reshape([20000.0, 60000.0, 100000.0], nλ, nθ, nd + 1)
    p_full = reshape([40000.0, 80000.0], nλ, nθ, nd)
    temperature = reshape([250.0, 285.0], nλ, nθ, nd)
    q_sat =
        [lscale_saturation_specific_humidity(atmo, temperature[k], p_full[k]) for k = 1:nd]
    humidity = reshape(1.15 .* q_sat, nλ, nθ, nd)
    initial_temperature = copy(temperature)
    initial_humidity = copy(humidity)
    temperature_tendency = zeros(nλ, nθ, nd)
    humidity_tendency = zeros(nλ, nθ, nd)
    liquid_water_content = zeros(nλ, nθ, nd)
    precipitation = zeros(nλ, nθ, 1)

    JGCM.Atmos_Param_Module.Lscale_Cond!(
        atmo,
        temperature,
        humidity,
        p_full,
        p_half,
        effective_dt,
        heating_scale,
        temperature_tendency,
        humidity_tendency,
        liquid_water_content,
        precipitation,
    )
    return (;
        atmo,
        p_half,
        initial_temperature,
        initial_humidity,
        temperature,
        humidity,
        temperature_tendency,
        humidity_tendency,
        liquid_water_content,
        precipitation,
    )
end

@testset "Large-scale condensation signs, units, and heating scale" begin
    effective_dt = 600.0
    full_heating = run_lscale_column(effective_dt, 1.0)

    @test all(full_heating.humidity_tendency .< 0.0)
    @test all(full_heating.temperature_tendency .> 0.0)
    @test full_heating.liquid_water_content ≈ -full_heating.humidity_tendency
    @test full_heating.precipitation[1] > 0.0

    layer_mass = diff(vec(full_heating.p_half)) ./ full_heating.atmo.grav
    expected_precipitation = -sum(vec(full_heating.humidity_tendency) .* layer_mass)
    @test full_heating.precipitation[1] ≈ expected_precipitation rtol = 1.0e-12

    energy_residual =
        full_heating.atmo.cp_air .* full_heating.temperature_tendency .+
        full_heating.atmo.Lv .* full_heating.humidity_tendency
    @test maximum(abs, energy_residual) < 1.0e-10

    weak_heating = run_lscale_column(effective_dt, 0.25)
    no_heating = run_lscale_column(effective_dt, 0.0)
    array_heating = run_lscale_column(effective_dt, fill(0.25, 1, 1))
    weak_energy_residual =
        weak_heating.atmo.cp_air .* weak_heating.temperature_tendency .+
        0.25 .* weak_heating.atmo.Lv .* weak_heating.humidity_tendency
    @test maximum(abs, weak_energy_residual) < 1.0e-10
    @test all(weak_heating.humidity_tendency .< full_heating.humidity_tendency)
    @test all(weak_heating.precipitation .> full_heating.precipitation)
    @test all(iszero, no_heating.temperature_tendency)
    @test all(no_heating.humidity_tendency .< weak_heating.humidity_tendency)
    @test all(no_heating.precipitation .> weak_heating.precipitation)
    @test array_heating.temperature_tendency ≈ weak_heating.temperature_tendency
    @test array_heating.humidity_tendency ≈ weak_heating.humidity_tendency
    @test array_heating.precipitation ≈ weak_heating.precipitation
end

@testset "Large-scale condensation finite working-state adjustment" begin
    short_step = run_lscale_column(600.0, 0.2)
    long_step = run_lscale_column(1200.0, 0.2)
    @test short_step.humidity_tendency .* 600.0 ≈ long_step.humidity_tendency .* 1200.0
    @test short_step.temperature_tendency .* 600.0 ≈
          long_step.temperature_tendency .* 1200.0
    @test short_step.precipitation .* 600.0 ≈ long_step.precipitation .* 1200.0
    @test short_step.temperature ≈
          short_step.initial_temperature .+ 600.0 .* short_step.temperature_tendency
    @test short_step.humidity ≈
          short_step.initial_humidity .+ 600.0 .* short_step.humidity_tendency
end

@testset "Large-scale condensation ignores negative spectral undershoots" begin
    nλ, nθ, nd = 1, 1, 2
    atmo = lscale_test_atmosphere(nλ, nθ, nd)
    p_half = reshape([20000.0, 60000.0, 100000.0], nλ, nθ, nd + 1)
    p_full = reshape([40000.0, 80000.0], nλ, nθ, nd)
    temperature = reshape([250.0, 285.0], nλ, nθ, nd)
    humidity = reshape([-1.0e-5, 0.0], nλ, nθ, nd)
    temperature_tendency = fill(NaN, nλ, nθ, nd)
    humidity_tendency = fill(NaN, nλ, nθ, nd)
    liquid_water_content = fill(NaN, nλ, nθ, nd)
    precipitation = zeros(nλ, nθ, 1)

    JGCM.Atmos_Param_Module.Lscale_Cond!(
        atmo,
        temperature,
        humidity,
        p_full,
        p_half,
        600.0,
        0.2,
        temperature_tendency,
        humidity_tendency,
        liquid_water_content,
        precipitation,
    )

    @test all(iszero, temperature_tendency)
    @test all(iszero, humidity_tendency)
    @test all(iszero, liquid_water_content)
    @test all(iszero, precipitation)
end

@testset "Shared Betts-Miller saturation derivative" begin
    thermo = JGCM.Atmos_Param_Module
    for temperature in (240.0, 263.16, 290.0)
        _, derivative = thermo._bm_saturation_vapor_pressure_and_derivative(temperature)
        step = 1.0e-3
        finite_difference =
            (
                Betts_Miller_Saturation_Vapor_Pressure(temperature + step) -
                Betts_Miller_Saturation_Vapor_Pressure(temperature - step)
            ) / (2.0 * step)
        @test derivative ≈ finite_difference rtol = 2.0e-7
    end
end
