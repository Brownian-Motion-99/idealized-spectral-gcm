using JGCM

function bm_test_atmosphere(nd)
    return Atmo_Data(
        "bm_test",
        1,
        1,
        nd,
        false,
        false,
        false,
        false,
        [0.0];
        radius = 6.371e6,
    )
end

@testset "Betts-Miller thermodynamics and validation" begin
    @test Betts_Miller_Saturation_Vapor_Pressure(240.0) > 0
    @test Betts_Miller_Saturation_Vapor_Pressure(263.16) > 0
    @test Betts_Miller_Saturation_Vapor_Pressure(290.0) > 0
    @test Betts_Miller_Saturation_Vapor_Pressure(290.0) >
          Betts_Miller_Saturation_Vapor_Pressure(263.16) >
          Betts_Miller_Saturation_Vapor_Pressure(240.0)
    @test isapprox(
        Betts_Miller_Saturation_Vapor_Pressure(253.16 - 1.0e-6),
        Betts_Miller_Saturation_Vapor_Pressure(253.16 + 1.0e-6);
        rtol = 1.0e-6,
    )
    @test isapprox(
        Betts_Miller_Saturation_Vapor_Pressure(273.16 - 1.0e-6),
        Betts_Miller_Saturation_Vapor_Pressure(273.16 + 1.0e-6);
        rtol = 1.0e-6,
    )

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
    invalid_humidity[3] = -1.0e-6
    @test_throws ArgumentError Betts_Miller_Column(
        state,
        atmo,
        temperature,
        invalid_humidity,
        p_full,
        p_half,
    )
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
    humidity = [5.0e-4, 2.0e-3, 5.0e-3, 1.0e-2, 1.6e-2, 2.0e-2]
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
    rs = epsilon * Betts_Miller_Saturation_Vapor_Pressure(300.0) / p_full[end]
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
        humidity[i, j, :] .= [5.0e-4, 2.0e-3, 5.0e-3, 1.0e-2, 1.6e-2, 2.0e-2]
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

    dyn_data = Dyn_Data("bm_storage", 1, 2, 4, 2, nd, 1)
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
    dyn = Dyn_Data("bm_interface", num_fourier, num_spherical, nλ, nθ, nd, 1)

    p_half_column = collect(range(0.0, 1.0e5; length = nd + 1))
    p_full_column = 0.5 .* (p_half_column[1:end-1] .+ p_half_column[2:end])
    temperature_column = [210.0, 225.0, 240.0, 255.0, 275.0, 300.0]
    humidity_column = [5.0e-4, 2.0e-3, 5.0e-3, 1.0e-2, 1.6e-2, 2.0e-2]
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

    JGCM.Atmos_Param_Module.Spectral_Physics!(config, mesh, vert, atmo, dyn, semi, params)
    @test dyn.grid_δt ≈ dyn.grid_bm_t_tendency + dyn.grid_lrf_tendency
    @test dyn.grid_δq ≈ dyn.grid_bm_q_tendency
    @test dyn.grid_precip ≈ dyn.grid_bm_precip
    @test dyn.grid_lrf_tendency ≈ 2.0 .* dyn.grid_q_c ./ config.day_to_sec
    @test !isapprox(dyn.grid_lrf_tendency, 2.0 .* dyn.grid_q_p ./ config.day_to_sec)

    params["do_Lscale_Cond"] = true
    @test_throws ErrorException JGCM.Atmos_Param_Module.Spectral_Physics!(
        config,
        mesh,
        vert,
        atmo,
        dyn,
        semi,
        params,
    )
    params["do_Lscale_Cond"] = false

    semi.integrator.init_step = false
    params["BM_state"] = Betts_Miller_State(nd; tau = 1000.0)
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
