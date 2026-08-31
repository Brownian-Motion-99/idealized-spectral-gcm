using Test
using JGCM

@testset "Moist Held-Suarez equilibrium profile" begin
    nλ, nθ, nd = 1, 3, 2
    sin_latitude = [-1.0, 0.0, 1.0]
    atmo = Atmo_Data(
        "hs_profile",
        nλ,
        nθ,
        nd,
        false,
        false,
        false,
        false,
        sin_latitude;
        radius=6371.0e3,
    )

    temperature = fill(250.0, nλ, nθ, nd)
    p_half = zeros(nλ, nθ, nd + 1)
    p_half[:, :, end] .= 1.0e5
    p_full = zeros(nλ, nθ, nd)
    p_full[:, :, 1] .= 1.0e3
    p_full[:, :, 2] .= 1.0e5
    equilibrium_temperature = similar(temperature)
    params = Dict{String,Any}(
        "σ_b" => 0.7,
        "k_a" => 1.0 / 40.0,
        "k_s" => 1.0 / 4.0,
        "k_f" => 1.0,
        "T_equator" => 294.0,
        "T_stratosphere" => 200.0,
        "ΔT_y" => 65.0,
        "Δθ_z" => 10.0,
    )

    JGCM.Atmos_Param_Module.Newtonian_Relaxation!(
        atmo,
        600,
        86400,
        sin_latitude,
        p_half,
        p_full,
        temperature,
        equilibrium_temperature,
        params,
    )

    @test equilibrium_temperature[1, 2, 2] ≈ 294.0
    @test equilibrium_temperature[1, 1, 2] ≈ 229.0
    @test equilibrium_temperature[1, 3, 2] ≈ 229.0
    @test all(equilibrium_temperature[:, :, 1] .== 200.0)

end

@testset "Held-Suarez uses the current leapfrog state" begin
    num_fourier, nθ, nd = 3, 12, 2
    num_spherical = num_fourier + 1
    nλ = 2 * nθ
    radius = 6371.0e3
    mesh = Spectral_Spherical_Mesh(
        num_fourier,
        num_spherical,
        nλ,
        nθ,
        nd,
        radius,
    )
    vert = Vert_Coordinate(
        nλ,
        nθ,
        nd,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )
    atmo = Atmo_Data(
        "hs_coupling",
        nλ,
        nθ,
        nd,
        false,
        false,
        false,
        false,
        mesh.sinθ;
        radius=radius,
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
    dyn = Dyn_Data(
        "hs_coupling",
        num_fourier,
        num_spherical,
        nλ,
        nθ,
        nd,
    )

    dyn.grid_ps_c .= 1.0e5
    dyn.grid_ps_p .= 1.0e5
    dyn.grid_p_half[:, :, 1] .= 0.0
    dyn.grid_p_half[:, :, 2] .= 5.0e4
    dyn.grid_p_half[:, :, 3] .= 1.0e5
    dyn.grid_p_full[:, :, 1] .= 2.5e4
    dyn.grid_p_full[:, :, 2] .= 7.5e4
    dyn.grid_u_c .= 12.0
    dyn.grid_v_c .= -6.0
    dyn.grid_t_c .= 300.0
    dyn.grid_ps_n .= dyn.grid_ps_c
    dyn.grid_u_n .= dyn.grid_u_c
    dyn.grid_v_n .= dyn.grid_v_c
    dyn.grid_t_n .= dyn.grid_t_c
    dyn.grid_u_p .= 120.0
    dyn.grid_v_p .= 60.0
    dyn.grid_t_p .= 240.0

    params = Dict{String,Any}(
        "do_Sensible_Heating" => false,
        "do_Surface_Evaporation" => false,
        "do_Implicit_PBL_Scheme" => false,
        "do_HS_Forcing" => true,
        "do_Betts_Miller" => false,
        "do_Lscale_Cond" => false,
        "do_LRF" => false,
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
        name="hs_coupling",
        model_type=:PrimitiveEquation,
        num_fourier=num_fourier,
        nθ=nθ,
        nd=nd,
        radius=radius,
        omega=7.292e-5,
        grav=9.8,
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
        moisture_processes=false,
        initial_condition=:Moist_Spinup,
        output_path="/tmp",
        output_filename="/tmp/hs_coupling.nc",
        logger="/tmp/hs_coupling.log",
        vars_to_output=Symbol[],
        output_interval=600,
        physics_params=params,
    )

    current_u = copy(dyn.grid_u_c)
    current_v = copy(dyn.grid_v_c)
    current_t = copy(dyn.grid_t_c)
    JGCM.Atmos_Param_Module.Spectral_Physics!(
        config,
        mesh,
        vert,
        atmo,
        dyn,
        semi,
        params,
    )

    sigma = 0.75
    damping_rate = (1.0 / 86400.0) * (sigma - 0.7) / (1.0 - 0.7)
    @test all((dyn.grid_u_n[:, :, 2] .- current_u[:, :, 2]) ./ 600.0 .≈
              -damping_rate * 12.0)
    @test all((dyn.grid_v_n[:, :, 2] .- current_v[:, :, 2]) ./ 600.0 .≈
              damping_rate * 6.0)
    @test !all((dyn.grid_u_n[:, :, 2] .- current_u[:, :, 2]) ./ 600.0 .≈
               -damping_rate * 120.0)

    equatorial_index = argmin(abs.(mesh.θc))
    sin_latitude = sin(mesh.θc[equatorial_index])
    cos_latitude = cos(mesh.θc[equatorial_index])
    p_norm = 0.75
    expected_equilibrium =
        (294.0 - 65.0 * sin_latitude^2 - 10.0 * cos_latitude^2 * log(p_norm)) *
        p_norm^atmo.kappa
    thermal_rate =
        (1.0 / 40.0 +
         (1.0 / 4.0 - 1.0 / 40.0) *
         (sigma - 0.7) / (1.0 - 0.7) * cos_latitude^4) /
        86400.0
    frictional_heating =
        -((12.0 - 0.5 * damping_rate * 12.0 * 600) * (-damping_rate * 12.0) +
          (-6.0 + 0.5 * damping_rate * 6.0 * 600) * (damping_rate * 6.0)) /
        atmo.cp_air
    post_friction_temperature = 300.0 + 600.0 * frictional_heating
    expected_dt =
        frictional_heating -
        thermal_rate * (post_friction_temperature - expected_equilibrium)
    @test dyn.grid_t_eq[1, equatorial_index, 2] ≈ expected_equilibrium
    @test (dyn.grid_t_n[1, equatorial_index, 2] - current_t[1, equatorial_index, 2]) /
          600.0 ≈ expected_dt
    @test dyn.grid_u_c == current_u
    @test dyn.grid_v_c == current_v
    @test dyn.grid_t_c == current_t

    workspace = params["Physics_workspace"]::Physics_Workspace
    kinetic_energy_loss =
        0.5 * (
            current_u[1, equatorial_index, 2]^2 +
            current_v[1, equatorial_index, 2]^2 -
            workspace.grid_u[1, equatorial_index, 2]^2 -
            workspace.grid_v[1, equatorial_index, 2]^2
        )
    friction_only_temperature = 300.0 + kinetic_energy_loss / atmo.cp_air
    expected_final_temperature =
        friction_only_temperature +
        600.0 * thermal_rate * (expected_equilibrium - friction_only_temperature)
    @test workspace.grid_t[1, equatorial_index, 2] ≈ expected_final_temperature
end
