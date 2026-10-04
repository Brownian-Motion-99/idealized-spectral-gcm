using JGCM

function bm_test_atmosphere(nd; use_virtual_temperature = false)
    return Atmo_Data(
        "bm_test",
        1,
        1,
        nd,
        false,
        false,
        false,
        use_virtual_temperature,
        [0.0];
        radius = 6.371e6,
    )
end

@testset "Betts-Miller dry ascent and LCL boundaries" begin
    nd = 6
    atmo = bm_test_atmosphere(nd; use_virtual_temperature = true)
    state = Betts_Miller_State(nd)
    epsilon = atmo.rdgas / atmo.rvgas
    p_half = [10000.0, 25000.0, 40000.0, 55000.0, 70000.0, 85000.0, 100000.0]
    p_full = 0.5 .* (p_half[1:end-1] .+ p_half[2:end])
    temperature = 300.0 .* (p_full ./ p_full[end]) .^ atmo.kappa

    # No in-domain condensation: preserve analytic dry ascent, even for r0=0.
    for q0 in (0.0, 1.0e-12), top_pressure in (0.0, 10000.0)
        interfaces = copy(p_half)
        interfaces[1] = top_pressure
        pressures = 0.5 .* (interfaces[1:end-1] .+ interfaces[2:end])
        dry_temperature = 300.0 .* (pressures ./ pressures[end]) .^ atmo.kappa
        humidity = fill(q0, nd)
        result = Betts_Miller_Column(state, atmo, dry_temperature, humidity, pressures, interfaces)
        @test !result.active
        @test result.lcl == result.lfc == result.lzb == 0
        @test result.cape == result.precipitation == 0.0
        @test abs(result.cin) < 1.0e-10
        @test result.parcel_temperature ≈ dry_temperature atol = 1.0e-12 rtol = 0
        @test all(==(q0 / (1 - q0)), result.parcel_mixing_ratio)
        @test all(iszero, result.temperature_tendency)
        @test all(iszero, result.humidity_tendency)
    end

    # This parcel condenses between levels 1 and 2 (pLCL ≈ 18540.835 Pa).
    # All levels below the LCL must stay dry, including level 2.
    q0 = 1.0e-6
    humidity = min.(q0, Saturation_Specific_Humidity.(temperature, p_full, epsilon))
    first_level = Betts_Miller_Column(state, atmo, temperature, humidity, p_full, p_half)
    @test first_level.lcl == 1
    @test first_level.parcel_temperature[2:end] ≈ temperature[2:end] atol = 1.0e-12 rtol = 0
    @test all(==(q0 / (1 - q0)), first_level.parcel_mixing_ratio[2:end])

    # Construct condensation exactly at a full level from dry-adiabatic
    # thermodynamics, independently of the LCL bisection routine.
    for k_lcl in (1, 3, nd)
        r0 = Saturation_Mixing_Ratio(temperature[k_lcl], p_full[k_lcl], epsilon)
        q0 = r0 / (1.0 + r0)
        humidity = min.(q0, Saturation_Specific_Humidity.(temperature, p_full, epsilon))
        humidity[end] = q0
        result = Betts_Miller_Column(state, atmo, temperature, humidity, p_full, p_half)
        @test result.lcl == k_lcl
        @test result.parcel_temperature[k_lcl:end] ≈ temperature[k_lcl:end] atol = 1.0e-12 rtol = 0
        @test result.parcel_mixing_ratio[k_lcl] ≈ r0 rtol = 1.0e-13
        @test all(isfinite, result.parcel_temperature)
        @test result.lfc == k_lcl - 1
    end

    # Below the LCL, equal environmental and parcel humidity makes this dry
    # adiabat neutrally buoyant. Substituting saturation humidity would give
    # spurious buoyancy and an incorrect signed contribution to CIN.
    q0 = 0.006
    humidity = min.(q0, Saturation_Specific_Humidity.(temperature, p_full, epsilon))
    conserved = Betts_Miller_Column(state, atmo, temperature, humidity, p_full, p_half)
    @test conserved.lcl == 4
    @test all(==(q0 / (1 - q0)), conserved.parcel_mixing_ratio[5:6])
    @test all(conserved.parcel_saturation_mixing_ratio[5:6] .> conserved.parcel_mixing_ratio[5:6])
    @test abs(conserved.cin) < 1.0e-10

    # A warm layer below the LCL gives an analytic CIN increment. Its virtual
    # correction must use conserved humidity, rather than saturation humidity.
    warm_layer = copy(temperature)
    warm_layer[5] += 2.0
    for virtual in (false, true)
        atmosphere = bm_test_atmosphere(nd; use_virtual_temperature = virtual)
        result = Betts_Miller_Column(state, atmosphere, warm_layer, humidity, p_full, p_half)
        coefficient = virtual ? atmo.rvgas / atmo.rdgas - 1 : 0.0
        expected_cin = atmo.rdgas * 2.0 * (1 + coefficient * q0) * log(p_half[6] / p_half[5])
        @test result.cin ≈ expected_cin atol = 1.0e-10 rtol = 0
    end

    # Two levels are sufficient for a valid LCL at the top full level.
    pressures = [70000.0, 95000.0]
    interfaces = [60000.0, 80000.0, 100000.0]
    dry_temperature = 300.0 .* (pressures ./ pressures[end]) .^ atmo.kappa
    r0 = Saturation_Mixing_Ratio(dry_temperature[1], pressures[1], epsilon)
    two_level = Betts_Miller_Column(
        Betts_Miller_State(2), bm_test_atmosphere(2; use_virtual_temperature = true),
        dry_temperature, fill(r0 / (1 + r0), 2), pressures, interfaces,
    )
    @test two_level.lcl == 1
    @test two_level.parcel_temperature ≈ dry_temperature atol = 1.0e-12 rtol = 0
    @test two_level.cape == 0.0
end

@testset "Betts-Miller virtual buoyancy and grid consistency" begin
    fixture = joinpath(@__DIR__, "fixtures", "betts_miller_virtual_column.tsv")
    rows = [parse.(Float64, split(line)) for line in readlines(fixture) if !startswith(line, "#")]
    data = reduce(hcat, rows)
    p_full, temperature, humidity = data[1, :], data[4, :], data[5, :]
    p_half = vcat(data[2, :], data[3, end])
    nd = length(p_full)
    state = Betts_Miller_State(nd)
    moist_atmo = bm_test_atmosphere(nd; use_virtual_temperature = true)
    moist = Betts_Miller_Column(state, moist_atmo, temperature, humidity, p_full, p_half)
    dry = Betts_Miller_Column(state, bm_test_atmosphere(nd), temperature, humidity, p_full, p_half)

    # Independent reference: unmodified Isca qe_moist_convection.F90, using
    # matched constants and double precision (details in fixtures/README.md).
    @test moist.cape ≈ 454.23657013397872 atol = 1.0e-8 rtol = 0
    @test moist.lcl == nd
    @test moist.lfc == nd - 1
    @test moist.lzb == 9
    @test dry.cape == 0.0
    # Positive CAPE alone does not satisfy the thermal/drying criteria.
    @test !moist.active
    @test moist.precipitation == 0.0

    # Exercise the threaded public entry point with both flag settings on a
    # column where virtual buoyancy changes the adjustment depth and rates.
    nd = 3
    state = Betts_Miller_State(nd)
    pressures = [75000.0, 85000.0, 95000.0]
    interfaces = [70000.0, 80000.0, 90000.0, 100000.0]
    temperatures = [293.0, 294.0, 300.0]
    humidities = [0.01, 0.0175, Saturation_Specific_Humidity(300.0, 95000.0, moist_atmo.rdgas / moist_atmo.rvgas)]
    grid(field) = repeat(reshape(field, 1, 1, length(field)), 2, 3, 1)
    results = []
    for virtual in (false, true)
        atmo = bm_test_atmosphere(nd; use_virtual_temperature = virtual)
        column = Betts_Miller_Column(state, atmo, temperatures, humidities, pressures, interfaces)
        @test column.active
        push!(results, column)
        dt, dq, rain = zeros(2, 3, nd), zeros(2, 3, nd), zeros(2, 3, 1)
        Betts_Miller!(state, atmo, grid(temperatures), grid(humidities), grid(pressures), grid(interfaces), dt, dq, rain)
        @test dt ≈ grid(column.temperature_tendency)
        @test dq ≈ grid(column.humidity_tendency)
        @test all(==(column.precipitation), rain)
    end
    @test results[1].lzb == 2
    @test results[2].lzb == 1
    @test results[1].cape != results[2].cape
    @test !isapprox(results[1].temperature_tendency, results[2].temperature_tendency)

    # A saturated parcel can remain buoyant all the way to a zero-pressure top.
    atmo = bm_test_atmosphere(2; use_virtual_temperature = true)
    epsilon = atmo.rdgas / atmo.rvgas
    result = Betts_Miller_Column(
        Betts_Miller_State(2), atmo, [260.0, 300.0],
        [0.001, Saturation_Specific_Humidity(300.0, 95000.0, epsilon)],
        [70000.0, 95000.0], [0.0, 80000.0, 100000.0],
    )
    @test result.lzb == result.lfc == 1
    @test isfinite(result.cape) && result.cape > 0
end

@testset "Betts-Miller thermodynamics and validation" begin
    @test isapprox(
        Saturation_Vapor_Pressure(253.16 - 1.0e-6),
        Saturation_Vapor_Pressure(253.16 + 1.0e-6);
        rtol = 1.0e-6,
    )
    @test isapprox(
        Saturation_Vapor_Pressure(273.16 - 1.0e-6),
        Saturation_Vapor_Pressure(273.16 + 1.0e-6);
        rtol = 1.0e-6,
    )

    temperature = 300.0
    pressure = 100_000.0
    epsilon = 287.0 / 461.5
    vapor_pressure = Saturation_Vapor_Pressure(temperature)
    expected_mixing_ratio = epsilon * vapor_pressure / (pressure - vapor_pressure)
    expected_specific_humidity =
        epsilon * vapor_pressure / (pressure - (1.0 - epsilon) * vapor_pressure)
    @test Saturation_Mixing_Ratio(temperature, pressure, epsilon) ≈
          expected_mixing_ratio rtol = 1.0e-14
    @test Saturation_Specific_Humidity(temperature, pressure, epsilon) ≈
          expected_specific_humidity rtol = 1.0e-14
    @test Saturation_Specific_Humidity(temperature, pressure, epsilon) ≈
          expected_mixing_ratio / (1.0 + expected_mixing_ratio) rtol = 1.0e-14
    @test_throws DomainError Saturation_Mixing_Ratio(temperature, vapor_pressure, epsilon)
    @test_throws DomainError Saturation_Specific_Humidity(temperature, vapor_pressure, epsilon)
    @test_throws ArgumentError Saturation_Mixing_Ratio(temperature, pressure, 1.0)

    @test_throws ArgumentError Betts_Miller_State(6; tau = 0.0)
    @test_throws ArgumentError Betts_Miller_State(6; relative_humidity = 0.0)

    nd = 6
    atmo = bm_test_atmosphere(nd)
    state = Betts_Miller_State(nd)
    p_half = [10000.0, 25000.0, 40000.0, 55000.0, 70000.0, 85000.0, 100000.0]
    p_full = 0.5 .* (p_half[1:end-1] .+ p_half[2:end])
    temperature = [230.0, 245.0, 260.0, 275.0, 290.0, 300.0]
    humidity = fill(1.0e-5, nd)

    @test_throws DimensionMismatch Betts_Miller_Column(
        state,
        atmo,
        temperature,
        humidity,
        p_full,
        p_half[1:end-1],
    )
    invalid_humidity = copy(humidity)
    invalid_humidity[3] = 1.0
    @test_throws ArgumentError Betts_Miller_Column(
        state,
        atmo,
        temperature,
        invalid_humidity,
        p_full,
        p_half,
    )
    invalid_humidity[3] = NaN
    @test_throws ArgumentError Betts_Miller_Column(
        state,
        atmo,
        temperature,
        invalid_humidity,
        p_full,
        p_half,
    )

    undershot_humidity = copy(humidity)
    undershot_humidity[3] = -1.3302935592021637e-5
    undershot = Betts_Miller_Column(
        state,
        atmo,
        temperature,
        undershot_humidity,
        p_full,
        p_half,
    )
    @test undershot_humidity[3] == -1.3302935592021637e-5
    @test all(undershot.reference_humidity .>= 0.0)
    @test all(isfinite, undershot.temperature_tendency)
    @test all(isfinite, undershot.humidity_tendency)

    @test_throws ArgumentError Betts_Miller_Column(
        state,
        atmo,
        temperature,
        humidity,
        reverse(p_full),
        p_half,
    )
end

@testset "Betts-Miller column physics" begin
    nd = 6
    atmo = bm_test_atmosphere(nd)
    state = Betts_Miller_State(nd; tau = 7200.0, relative_humidity = 0.8)
    p_half = [10000.0, 25000.0, 40000.0, 55000.0, 70000.0, 85000.0, 100000.0]
    p_full = 0.5 .* (p_half[1:end-1] .+ p_half[2:end])

    stable = Betts_Miller_Column(
        state,
        atmo,
        [230.0, 245.0, 260.0, 275.0, 290.0, 300.0],
        fill(1.0e-5, nd),
        p_full,
        p_half,
    )
    @test !stable.active
    @test stable.cape == 0.0
    @test all(iszero, stable.temperature_tendency)
    @test all(iszero, stable.humidity_tendency)
    @test stable.precipitation == 0.0

    temperature = [210.0, 225.0, 240.0, 255.0, 275.0, 300.0]
    humidity = [6.0e-4, 2.4e-3, 6.0e-3, 1.2e-2, 1.92e-2, 2.4e-2]
    convective = Betts_Miller_Column(state, atmo, temperature, humidity, p_full, p_half)

    @test convective.active
    @test convective.cape > 0
    @test 1 <= convective.lzb <= convective.lfc <= nd
    @test convective.precipitation > 0
    @test any(convective.temperature_tendency .> 0)
    @test sum(convective.humidity_tendency) < 0
    @test all(iszero, convective.temperature_tendency[1:convective.lzb-1])
    @test all(iszero, convective.humidity_tendency[1:convective.lzb-1])

    layer_mass = diff(p_half) ./ atmo.grav
    energy_residual = sum(
        (
            atmo.cp_air .* convective.temperature_tendency .+
            atmo.Lv .* convective.humidity_tendency
        ) .* layer_mass,
    )
    moisture_precipitation = -sum(convective.humidity_tendency .* layer_mass)
    @test abs(energy_residual) < 1.0e-8
    @test moisture_precipitation ≈ convective.precipitation rtol = 1.0e-12

    # A saturated surface parcel must use the explicit surface saturation branch.
    epsilon = atmo.rdgas / atmo.rvgas
    rs = Saturation_Mixing_Ratio(300.0, p_full[end], epsilon)
    saturated_humidity = copy(humidity)
    saturated_humidity[end] = rs / (1.0 + rs)
    saturated =
        Betts_Miller_Column(state, atmo, temperature, saturated_humidity, p_full, p_half)
    @test saturated.lcl == nd
end

@testset "Betts-Miller grid buffers and output registration" begin
    nλ, nθ, nd = 2, 2, 6
    atmo = Atmo_Data(
        "bm_grid_test",
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
    state = Betts_Miller_State(nd)
    p_half_column = [10000.0, 25000.0, 40000.0, 55000.0, 70000.0, 85000.0, 100000.0]
    p_full_column = 0.5 .* (p_half_column[1:end-1] .+ p_half_column[2:end])
    temperature = zeros(nλ, nθ, nd)
    humidity = zeros(nλ, nθ, nd)
    p_full = zeros(nλ, nθ, nd)
    p_half = zeros(nλ, nθ, nd + 1)
    for j = 1:nθ, i = 1:nλ
        temperature[i, j, :] .= [210.0, 225.0, 240.0, 255.0, 275.0, 300.0]
        humidity[i, j, :] .= [6.0e-4, 2.4e-3, 6.0e-3, 1.2e-2, 1.92e-2, 2.4e-2]
        p_full[i, j, :] .= p_full_column
        p_half[i, j, :] .= p_half_column
    end
    bm_dt = fill(NaN, nλ, nθ, nd)
    bm_dq = fill(NaN, nλ, nθ, nd)
    bm_precip = fill(NaN, nλ, nθ, 1)
    Betts_Miller!(
        state,
        atmo,
        temperature,
        humidity,
        p_full,
        p_half,
        bm_dt,
        bm_dq,
        bm_precip,
    )
    @test all(isfinite, bm_dt)
    @test all(isfinite, bm_dq)
    @test all(bm_precip .> 0)

    dyn_data = Dyn_Data("bm_storage", 1, 2, 4, 2, nd)
    live_data =
        JGCM.Variable_Mappings_Module.Get_Dyn_Var_Map(dyn_data, Val(:PrimitiveEquation))
    @test live_data[:bm_dt] === dyn_data.grid_bm_t_tendency
    @test live_data[:bm_dq] === dyn_data.grid_bm_q_tendency
    @test live_data[:bm_precip] === dyn_data.grid_bm_precip
    metadata = JGCM.Output_Mappings_Module.Get_Var_Info(Val(:PrimitiveEquation))
    @test metadata[:bm_dt].units == "K s-1"
    @test metadata[:bm_dq].units == "s-1"
    @test metadata[:bm_precip].units == "kg m-2 s-1"
end

@testset "Betts-Miller and LRF physics-interface accumulation" begin
    num_fourier, nθ, nd = 1, 16, 6
    num_spherical = num_fourier + 1
    nλ = 2nθ
    radius = 6.371e6
    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nλ, nθ, nd, radius)
    vert = Vert_Coordinate(
        nλ,
        nθ,
        nd,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )
    atmo = Atmo_Data(
        "bm_interface",
        nλ,
        nθ,
        nd,
        false,
        false,
        false,
        false,
        mesh.sinθ;
        radius = radius,
    )
    integrator = JGCM.Time_Integrator_Module.Filtered_Leapfrog(
        0.04,
        4,
        1.0e-4,
        mesh.laplacian_eig,
        0.5,
        600,
        true,
        0,
        1200,
    )
    semi = JGCM.Semi_Implicit_Module.Semi_Implicit_Solver(
        vert,
        atmo,
        integrator,
        1.0e5,
        fill(300.0, nd),
        mesh.wave_numbers,
    )
    dyn = Dyn_Data("bm_interface", num_fourier, num_spherical, nλ, nθ, nd)

    p_half_column = collect(range(0.0, 1.0e5; length = nd + 1))
    p_full_column = 0.5 .* (p_half_column[1:end-1] .+ p_half_column[2:end])
    temperature_column = [210.0, 225.0, 240.0, 255.0, 275.0, 300.0]
    humidity_column = [6.0e-4, 2.4e-3, 6.0e-3, 1.2e-2, 1.92e-2, 2.4e-2]
    for j = 1:nθ, i = 1:nλ
        dyn.grid_ps_c[i, j, 1] = 1.0e5
        dyn.grid_ps_p[i, j, 1] = 1.0e5
        dyn.grid_p_half[i, j, :] .= p_half_column
        dyn.grid_p_full[i, j, :] .= p_full_column
        dyn.grid_t_c[i, j, :] .= temperature_column
        dyn.grid_t_p[i, j, :] .= temperature_column
        dyn.grid_q_c[i, j, :] .= humidity_column
        dyn.grid_q_p[i, j, :] .= 0.5 .* humidity_column
    end

    lrf_matrix = zeros(nd, nd, nθ)
    for j = 1:nθ, k = 1:nd
        lrf_matrix[k, k, j] = 2.0
    end
    lrf_state = LRF_State(lrf_matrix, zeros(nλ, nθ, nd))
    bm_state = Betts_Miller_State(nd; tau = 7200.0)
    params = Dict{String,Any}(
        "do_Lscale_Cond" => false,
        "do_Sensible_Heating" => false,
        "do_Surface_Evaporation" => false,
        "do_Implicit_PBL_Scheme" => false,
        "do_HS_Forcing" => false,
        "do_Betts_Miller" => true,
        "BM_state" => bm_state,
        "condensation_heating_fraction" => 0.2,
        "do_LRF" => true,
        "LRF_state" => lrf_state,
    )
    config = Model_Config(
        name = "bm_interface",
        model_type = :PrimitiveEquation,
        num_fourier = num_fourier,
        nθ = nθ,
        nd = nd,
        radius = radius,
        omega = 7.292e-5,
        grav = 9.8,
        vert_coord_option = "even_sigma",
        vert_difference_option = "simmons_and_burridge",
        vert_ref_level_option = "second_centered_wts",
        Δt = 600,
        end_time = 1200,
        day_to_sec = 86400,
        damping_order = 4,
        damping_coef = 1.0e-4,
        robert_coef = 0.04,
        implicit_coef = 0.5,
        moisture_processes = true,
        initial_condition = :Moist_Spinup,
        output_path = "/tmp",
        output_filename = "/tmp/bm_interface.nc",
        logger = "/tmp/bm_interface.log",
        vars_to_output = Symbol[],
        output_interval = 600,
        physics_params = params,
    )

    current_t = copy(dyn.grid_t_c)
    current_q = copy(dyn.grid_q_c)
    current_u = copy(dyn.grid_u_c)
    current_v = copy(dyn.grid_v_c)
    current_spe_t = copy(dyn.spe_t_c)
    dyn.grid_ps_n .= dyn.grid_ps_c
    dyn.grid_u_n .= dyn.grid_u_c
    dyn.grid_v_n .= dyn.grid_v_c
    dyn.grid_t_n .= dyn.grid_t_c
    dyn.grid_q_n .= dyn.grid_q_c

    JGCM.Atmos_Param_Module.Spectral_Physics!(config, mesh, vert, atmo, dyn, semi, params)
    workspace = params["Physics_workspace"]::Physics_Workspace
    @test (dyn.grid_t_n .- current_t) ./ 600.0 ≈
          dyn.grid_bm_t_tendency + dyn.grid_lrf_tendency
    q_after_physics = current_q .+ 600.0 .* dyn.grid_bm_q_tendency
    @test dyn.grid_q_n .* dyn.grid_Δp ≈
          q_after_physics .* workspace.grid_Δp_before
    @test dyn.grid_precip ≈ dyn.grid_bm_precip
    @test dyn.grid_lrf_tendency ≈ 2.0 .* q_after_physics ./ config.day_to_sec
    @test !isapprox(dyn.grid_lrf_tendency, 2.0 .* dyn.grid_q_p ./ config.day_to_sec)
    @test dyn.grid_t_n ≈ workspace.grid_t
    @test dyn.grid_q_n ≈ workspace.grid_q
    @test dyn.grid_t_c == current_t
    @test dyn.grid_q_c == current_q
    @test dyn.grid_u_c == current_u
    @test dyn.grid_v_c == current_v
    @test dyn.spe_t_c == current_spe_t

    params["do_Lscale_Cond"] = true
    dyn.grid_u_n .= current_u
    dyn.grid_v_n .= current_v
    dyn.grid_t_n .= current_t
    dyn.grid_q_n .= current_q
    JGCM.Atmos_Param_Module.Spectral_Physics!(config, mesh, vert, atmo, dyn, semi, params)
    @test all(isfinite, dyn.grid_t_n)
    @test all(isfinite, dyn.grid_q_n)
    @test all(dyn.grid_precip .>= dyn.grid_bm_precip)
    @test (dyn.grid_t_n .- current_t) ./ 600.0 ≈
          dyn.grid_bm_t_tendency + workspace.grid_lscale_t_tendency +
          dyn.grid_lrf_tendency
    q_after_physics = current_q .+ 600.0 .* (
        dyn.grid_bm_q_tendency + workspace.grid_lscale_q_tendency
    )
    @test dyn.grid_q_n .* dyn.grid_Δp ≈
          q_after_physics .* workspace.grid_Δp_before
    @test dyn.grid_lrf_tendency ≈ 2.0 .* q_after_physics ./ config.day_to_sec
    @test dyn.grid_t_n ≈ workspace.grid_t
    @test dyn.grid_q_n ≈ workspace.grid_q
    @test dyn.grid_t_c == current_t
    @test dyn.grid_q_c == current_q

    params["do_Lscale_Cond"] = false

    semi.integrator.init_step = false
    params["BM_state"] = Betts_Miller_State(nd; tau = 1000.0)
    dyn.grid_t_n .= current_t
    dyn.grid_q_n .= current_q
    JGCM.Atmos_Param_Module.Spectral_Physics!(config, mesh, vert, atmo, dyn, semi, params)
    @test all(isfinite, dyn.grid_t_n)
    @test all(isfinite, dyn.grid_q_n)

    params["BM_state"] = Betts_Miller_State(nd; tau = 500.0)
    @test_throws ArgumentError JGCM.Atmos_Param_Module.Spectral_Physics!(
        config,
        mesh,
        vert,
        atmo,
        dyn,
        semi,
        params,
    )
end
