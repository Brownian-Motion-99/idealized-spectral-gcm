function bm_coupling_frame(dt, initial_step)
    num_fourier, nlat, nd = 1, 16, 30
    nlon, radius = 2nlat, 6.371e6
    mesh = Spectral_Spherical_Mesh(num_fourier, 2, nlon, nlat, nd, radius)
    vert = Vert_Coordinate(nlon, nlat, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts")
    # Finite-top hybrid coefficients reproduce the fixture interfaces at 1 bar.
    vert.ak .= 1000.0 .* (1 .- vert.bk)
    vert.Δak .= diff(vert.ak)
    vert.zero_top = false
    atmo = Atmo_Data("bm_coupling", nlon, nlat, nd, false, false, false, true, mesh.sinθ; radius)
    integrator = JGCM.Time_Integrator_Module.Filtered_Leapfrog(
        0.04, 4, 1.0e-4, mesh.laplacian_eig, 0.5, dt, initial_step, 0, 2dt,
    )
    semi = JGCM.Semi_Implicit_Module.Semi_Implicit_Solver(
        vert, atmo, integrator, 1.0e5, fill(300.0, nd), mesh.wave_numbers,
    )
    dyn = Dyn_Data("bm_coupling", num_fourier, 2, nlon, nlat, nd)
    dyn.grid_ps_n .= 1.0e5
    JGCM.Press_And_Geopot_Module.Pressure_Variables!(
        vert, dyn.grid_ps_n, dyn.grid_p_half, dyn.grid_Δp, dyn.grid_lnp_half,
        dyn.grid_p_full, dyn.grid_lnp_full,
    )
    deep, shallow = BMValidation.fixture("deep"), BMValidation.fixture("shallow")
    for j = 1:nlat, i = 1:nlon
        case = isodd(i) ? deep : shallow
        dyn.grid_t_n[i, j, :] .= case.temperature
        dyn.grid_q_n[i, j, :] .= case.humidity
        # Ensure condensation is present, including a layer modified by BM.
        dyn.grid_q_n[i, j, 2] = 1.05 * Saturation_Specific_Humidity(
            case.temperature[2], dyn.grid_p_full[i, j, 2], atmo.rdgas / atmo.rvgas,
        )
        if isodd(i)
            dyn.grid_q_n[i, j, 25] = 1.01 * Saturation_Specific_Humidity(
                case.temperature[25], dyn.grid_p_full[i, j, 25], atmo.rdgas / atmo.rvgas,
            )
        end
    end
    # Deliberately different current/previous states catch coupling to stale data.
    dyn.grid_t_c .= dyn.grid_t_n .+ 5.0
    dyn.grid_q_c .= 0.5 .* dyn.grid_q_n
    dyn.grid_q_p .= 0.25 .* dyn.grid_q_n
    params = Dict{String,Any}(
        "do_Betts_Miller" => true, "BM_state" => Betts_Miller_State(nd),
        "do_Lscale_Cond" => false, "condensation_heating_fraction" => 1.0,
    )
    config = Model_Config(
        name = "bm_coupling", model_type = :PrimitiveEquation, num_fourier = num_fourier,
        nθ = nlat, nd = nd, radius = radius, omega = 7.292e-5, grav = 9.8,
        vert_coord_option = "even_sigma", vert_difference_option = "simmons_and_burridge",
        vert_ref_level_option = "second_centered_wts", Δt = dt, end_time = 2dt,
        day_to_sec = 86400, damping_order = 4, damping_coef = 1.0e-4, robert_coef = 0.04,
        implicit_coef = 0.5, moisture_processes = true, initial_condition = :Moist_Spinup,
        output_path = "/tmp", output_filename = "/tmp/bm_coupling.nc", logger = "/tmp/bm_coupling.log",
        vars_to_output = Symbol[], output_interval = dt, physics_params = params,
    )
    return (; config, mesh, vert, atmo, dyn, semi, params)
end

@testset "Betts-Miller ordered coupling and bounded substeps" begin
    for dt in (600, 7200), initial_step in (true, false), condensation in (false, true)
        @testset "dt=$dt startup=$initial_step condensation=$condensation" begin
            frame = bm_coupling_frame(dt, initial_step)
            config, mesh, vert, atmo, dyn, semi, params = frame
            params["do_Lscale_Cond"] = condensation
            nsteps = initial_step ? 1 : 2
            total_dt = nsteps * dt
            initial_t, initial_q = copy(dyn.grid_t_n), copy(dyn.grid_q_n)
            initial_mass = copy(dyn.grid_Δp) ./ atmo.grav
            current_t, current_q, previous_q = copy(dyn.grid_t_c), copy(dyn.grid_q_c), copy(dyn.grid_q_p)
            initial_water = sum(initial_q .* initial_mass; dims = 3)
            initial_dry = sum((1 .- initial_q) .* initial_mass; dims = 3)
            first_results = [Betts_Miller_Column(params["BM_state"], atmo,
                initial_t[i, 1, :], initial_q[i, 1, :], dyn.grid_p_full[i, 1, :], dyn.grid_p_half[i, 1, :]) for i in 1:2]
            @test first_results[1].regime == :deep
            @test first_results[2].regime == :shallow

            # Independent sequence through the public column routine and
            # condensation routine; do not call the combined physics wrapper.
            expected_t, expected_q = copy(initial_t), copy(initial_q)
            expected_ps = copy(dyn.grid_ps_n)
            pf, ph = copy(dyn.grid_p_full), copy(dyn.grid_p_half)
            dp, logph, logpf = copy(dyn.grid_Δp), similar(ph), similar(pf)
            average_t, average_q = zeros(size(initial_t)), zeros(size(initial_q))
            average_bm_rain, average_rain = zeros(size(expected_ps)), zeros(size(expected_ps))
            ldt, ldq, liquid, rain = zeros(size(initial_t)), zeros(size(initial_q)), zeros(size(initial_q)), zeros(size(expected_ps))
            counterfactual_difference = 0.0
            for step = 1:nsteps
                q_before = copy(expected_q)
                dp_before = copy(dp)
                pre_bm_t, pre_bm_q = copy(expected_t), copy(expected_q)
                fill!(rain, 0.0)
                for j in axes(expected_t, 2), i in axes(expected_t, 1)
                    column = Betts_Miller_Column(params["BM_state"], atmo,
                        expected_t[i, j, :], expected_q[i, j, :], pf[i, j, :], ph[i, j, :])
                    expected_t[i, j, :] .+= dt .* column.temperature_tendency
                    expected_q[i, j, :] .+= dt .* column.humidity_tendency
                    average_t[i, j, :] .+= column.temperature_tendency ./ nsteps
                    average_q[i, j, :] .+= column.humidity_tendency ./ nsteps
                    average_bm_rain[i, j, 1] += column.precipitation / nsteps
                    rain[i, j, 1] = column.precipitation
                end
                if condensation
                    JGCM.Atmos_Param_Module.Lscale_Cond!(atmo, expected_t, expected_q, pf, ph,
                        dt, 1.0, ldt, ldq, liquid, rain)
                    if step == 1
                        counter_rain = zeros(size(rain))
                        counter_ldt, counter_ldq = similar(ldt), similar(ldq)
                        JGCM.Atmos_Param_Module.Lscale_Cond!(atmo, pre_bm_t, pre_bm_q, pf, ph,
                            dt, 1.0, counter_ldt, counter_ldq, similar(liquid), counter_rain)
                        counterfactual_difference = maximum(abs.(counter_ldq - ldq))
                    end
                end
                average_rain .+= rain ./ nsteps
                JGCM.Atmos_Param_Module.Dry_Air_Adjustment!(vert, expected_ps, q_before, expected_q, dp_before)
                JGCM.Press_And_Geopot_Module.Pressure_Variables!(vert, expected_ps, ph, dp, logph, pf, logpf)
            end
            JGCM.Atmos_Param_Module.Spectral_Physics!(config, mesh, vert, atmo, dyn, semi, params)
            @test dyn.grid_t_n ≈ expected_t rtol = 1.0e-14
            @test dyn.grid_q_n ≈ expected_q rtol = 1.0e-14
            @test dyn.grid_ps_n ≈ expected_ps rtol = 1.0e-14
            @test dyn.grid_bm_t_tendency ≈ average_t rtol = 1.0e-13
            @test dyn.grid_bm_q_tendency ≈ average_q rtol = 1.0e-13
            @test dyn.grid_bm_precip ≈ average_bm_rain rtol = 1.0e-13
            @test dyn.grid_precip ≈ average_rain rtol = 1.0e-13
            @test all(iszero, dyn.grid_bm_precip[2:2:end, :, :])
            @test any(!iszero, dyn.grid_bm_q_tendency[2:2:end, :, :])
            @test dyn.grid_t_c == current_t && dyn.grid_q_c == current_q && dyn.grid_q_p == previous_q
            if condensation
                @test counterfactual_difference > 1.0e-12
                @test any(dyn.grid_precip .> dyn.grid_bm_precip)
            else
                @test (dyn.grid_t_n - initial_t) ./ total_dt ≈ dyn.grid_bm_t_tendency rtol = 1.0e-12
                @test dyn.grid_precip == dyn.grid_bm_precip
            end
            # The final humidity budget must use final layer masses, because
            # precipitation changes pressure and the dry-air adjustment rescales q.
            final_mass = dyn.grid_Δp ./ atmo.grav
            final_water = sum(dyn.grid_q_n .* final_mass; dims = 3)
            final_dry = sum((1 .- dyn.grid_q_n) .* final_mass; dims = 3)
            @test final_dry ≈ initial_dry rtol = 2.0e-14
            @test final_water .+ total_dt .* dyn.grid_precip ≈ initial_water rtol = 2.0e-14
            @test all(isfinite, dyn.grid_q_n) && all(dyn.grid_q_n .>= 0)
        end
    end
end
