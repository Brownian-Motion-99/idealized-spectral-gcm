using Test
using JGCM

@testset "Post-dynamics moist-physics coupling" begin
    num_fourier, nθ, nd = 5, 12, 8
    num_spherical = num_fourier + 1
    nλ = 2nθ
    radius = 6.371e6
    mesh = Spectral_Spherical_Mesh(
        num_fourier, num_spherical, nλ, nθ, nd, radius,
    )
    vert = Vert_Coordinate(
        nλ, nθ, nd, "even_sigma", "simmons_and_burridge",
        "second_centered_wts",
    )
    physics_params = Dict{String,Any}(
        "do_mass_correction" => true,
        "do_energy_correction" => true,
        "do_water_correction" => true,
        "use_virtual_temperature" => true,
        "initial_humidity_floor" => 0.0,
        "do_Betts_Miller" => false,
        "do_Lscale_Cond" => true,
        "condensation_heating_fraction" => 1.0,
        "do_LRF" => false,
        "do_Sensible_Heating" => true,
        "C_H" => 0.0044,
        "do_Surface_Evaporation" => true,
        "C_E" => 0.0044,
        "do_Implicit_PBL_Scheme" => true,
        "C_D" => 0.0044,
        "PBL_Top_Mode" => :PressureLevel,
        "PBL_Top_Value" => 85000.0,
        "do_HS_Forcing" => true,
        "σ_b" => 0.7,
        "k_a" => 1.0 / 40.0,
        "k_s" => 1.0 / 4.0,
        "k_f" => 1.0,
        "T_equator" => 294.0,
        "T_stratosphere" => 200.0,
        "ΔT_y" => 65.0,
        "Δθ_z" => 10.0,
    )
    config = Model_Config(
        name = "coupling_smoke",
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
        end_time = 18_000,
        day_to_sec = 86_400,
        damping_order = 4,
        damping_coef = 1.15741e-4,
        robert_coef = 0.04,
        implicit_coef = 0.5,
        moisture_processes = true,
        num_tracers = 0,
        initial_condition = :Moist_Spinup,
        output_path = "/tmp",
        output_filename = "/tmp/coupling_smoke.nc",
        logger = "/tmp/coupling_smoke.log",
        vars_to_output = Symbol[],
        output_interval = 600,
        physics_params = physics_params,
    )
    atmo = Atmo_Data(
        "coupling_smoke", nλ, nθ, nd, true, true, true, true, mesh.sinθ;
        radius = radius, omega = config.omega, grav = config.grav,
    )
    dyn = Dyn_Data(
        "coupling_smoke", num_fourier, num_spherical, nλ, nθ, nd, 0,
    )
    Initialize_Atmos_State!(mesh, atmo, dyn, vert, config)
    integrator = JGCM.Time_Integrator_Module.Filtered_Leapfrog(
        config.robert_coef, config.damping_order, config.damping_coef,
        mesh.laplacian_eig, config.implicit_coef, config.Δt, true, 0,
        config.end_time,
    )
    semi = JGCM.Semi_Implicit_Module.Semi_Implicit_Solver(
        vert, atmo, integrator, 1.0e5, fill(300.0, nd), mesh.wave_numbers,
    )

    for step in 1:30
        JGCM.Driver.Step_Dynamics!(
            config, mesh, atmo, dyn, integrator, semi, vert, physics_params,
        )
        if integrator.init_step
            JGCM.Semi_Implicit_Module.Update_Init_Step!(semi)
        end
        integrator.time += config.Δt

        @test all(isfinite, dyn.grid_u_c)
        @test all(isfinite, dyn.grid_v_c)
        @test all(isfinite, dyn.grid_t_c)
        @test all(isfinite, dyn.grid_q_c)
        @test minimum(dyn.grid_t_c) > 0.0
        @test minimum(dyn.grid_q_c) >= -1.0e-14
    end

    q_from_spectral = similar(dyn.grid_q_c)
    Trans_Spherical_To_Grid!(mesh, dyn.spe_q_c, q_from_spectral)
    @test q_from_spectral ≈ dyn.grid_q_c rtol = 1.0e-12 atol = 1.0e-14

    t_from_spectral = similar(dyn.grid_t_c)
    Trans_Spherical_To_Grid!(mesh, dyn.spe_t_c, t_from_spectral)
    @test t_from_spectral ≈ dyn.grid_t_c rtol = 1.0e-12 atol = 1.0e-12

    lnps_from_spectral = similar(dyn.grid_ps_c)
    Trans_Spherical_To_Grid!(mesh, dyn.spe_lnps_c, lnps_from_spectral)
    @test exp.(lnps_from_spectral) ≈ dyn.grid_ps_c rtol = 1.0e-12 atol = 1.0e-9

    u_from_spectral = similar(dyn.grid_u_c)
    v_from_spectral = similar(dyn.grid_v_c)
    UV_Grid_From_Vor_Div!(
        mesh, dyn.spe_vor_c, dyn.spe_div_c, u_from_spectral, v_from_spectral,
    )
    @test u_from_spectral ≈ dyn.grid_u_c rtol = 1.0e-11 atol = 1.0e-11
    @test v_from_spectral ≈ dyn.grid_v_c rtol = 1.0e-11 atol = 1.0e-11
end
