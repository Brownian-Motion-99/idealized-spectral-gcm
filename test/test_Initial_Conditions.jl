using Test
using JGCM

@testset "Moist aquaplanet initial condition" begin
    num_fourier = 15
    num_spherical = num_fourier + 1
    nθ = 24
    nλ = 2 * nθ
    nd = 10
    radius = 6371.0e3

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
        "moist_initial_condition",
        nλ,
        nθ,
        nd,
        true,
        true,
        false,
        true,
        mesh.sinθ;
        radius=radius,
        omega=7.292e-5,
        grav=9.80,
    )
    dyn = Dyn_Data(
        "moist_initial_condition",
        num_fourier,
        num_spherical,
        nλ,
        nθ,
        nd,
    )
    config = Model_Config(
        name="moist_initial_condition",
        model_type=:PrimitiveEquation,
        num_fourier=num_fourier,
        nθ=nθ,
        nd=nd,
        radius=radius,
        omega=7.292e-5,
        grav=9.80,
        vert_coord_option="even_sigma",
        vert_difference_option="simmons_and_burridge",
        vert_ref_level_option="second_centered_wts",
        Δt=600,
        end_time=1200,
        day_to_sec=86400,
        damping_order=4,
        damping_coef=1.0e-4,
        robert_coef=0.04,
        implicit_coef=0.5,
        moisture_processes=true,
        initial_condition=:Moist_Spinup,
        output_path="/tmp",
        output_filename="/tmp/moist_initial_condition.nc",
        logger="/tmp/moist_initial_condition.log",
        vars_to_output=Symbol[],
        output_interval=600,
        physics_params=Dict{String,Any}(
            "do_mass_correction" => true,
            "do_energy_correction" => true,
            "do_water_correction" => false,
            "use_virtual_temperature" => true,
            "initial_humidity_floor" => 0.0,
        ),
    )

    equator_surface =
        JGCM.Initial_Conditions.Ullrich_Shallow_Basic_State(atmo, 0.0, 1.0e5)
    pole_surface =
        JGCM.Initial_Conditions.Ullrich_Shallow_Basic_State(atmo, pi / 2.0, 1.0e5)
    @test equator_surface.height ≈ 0.0 atol = 1.0e-10
    @test equator_surface.virtual_temperature ≈ 310.0 atol = 1.0e-10
    @test equator_surface.zonal_wind ≈ 0.0 atol = 1.0e-10
    @test pole_surface.virtual_temperature ≈ 240.0 atol = 1.0e-10

    Initialize_Atmos_State!(mesh, atmo, dyn, vert, config)

    @test all(isfinite, dyn.grid_t_c)
    @test all(isfinite, dyn.grid_q_c)
    @test all(isfinite, dyn.grid_u_c)
    @test all(isfinite, dyn.grid_v_c)
    @test dyn.grid_ps_c ≈ fill(1.0e5, size(dyn.grid_ps_c)) rtol = 1.0e-12
    @test minimum(dyn.grid_t_c) > 150.0
    @test maximum(dyn.grid_t_c) < 320.0
    @test maximum(dyn.grid_q_c) > 0.015
    @test minimum(dyn.grid_q_c) > -1.0e-12
    @test maximum(abs, dyn.grid_u_c) > 20.0
    @test maximum(abs, dyn.grid_v_c) < 1.1

    epsilon = atmo.rdgas / atmo.rvgas
    relative_humidity = similar(dyn.grid_q_c)
    for index in eachindex(relative_humidity)
        saturation_specific_humidity = Saturation_Specific_Humidity(
            dyn.grid_t_c[index],
            dyn.grid_p_full[index],
            epsilon,
        )
        relative_humidity[index] = dyn.grid_q_c[index] / saturation_specific_humidity
    end
    @test maximum(relative_humidity) < 1.0

    top_levels = findall(dyn.grid_p_full[1, 1, :] .< 1.0e4)
    @test !isempty(top_levels)
    @test maximum(abs, dyn.grid_q_c[:, :, top_levels]) < 1.0e-12

    @test dyn.spe_vor_p == dyn.spe_vor_c
    @test dyn.spe_div_p == dyn.spe_div_c
    @test dyn.spe_t_p == dyn.spe_t_c
    @test dyn.grid_t_p == dyn.grid_t_c
    @test dyn.grid_q_p == dyn.grid_q_c
    @test !hasproperty(dyn, :spe_q_c)
end
