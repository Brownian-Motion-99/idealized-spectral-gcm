using Test
using NCDatasets
using JGCM

include(joinpath(@__DIR__, "..", "post_processing", "Interpolator.jl"))

@testset "Output validation and final partial interval" begin
    mesh = Spectral_Spherical_Mesh(1, 2, 8, 4, 1, 6.371e6)
    atmo = Atmo_Data(
        "output", 8, 4, 1, true, true, false, false, mesh.sinθ;
        radius = mesh.radius,
    )
    dyn = Dyn_Data("output", 1, 2, 8, 4, 1)
    dyn.grid_u_c .= 10.0

    mktempdir() do dir
        filename = joinpath(dir, "barotropic.nc")
        requested = [:u]
        manager = Output_Manager(
            mesh, nothing, atmo, 0, 600, requested;
            filename,
            model_mode = :Barotropic,
            output_interval = 3600,
        )
        Update_Output!(manager, dyn, 600)
        dyn.grid_u_c .= 15.0
        Update_Output!(manager, dyn, 1200)
        Finalize_Output!(manager)

        @test requested == [:u]
        NCDataset(joinpath(dir, "barotropic_t0.nc"), "r") do ds
            @test ds["time"].var[:] ≈ [1200 / 86_400]
            @test all(ds["ua"].var[:, :, 1] .== 12.5)
        end
    end

    mktempdir() do dir
        @test_throws ArgumentError Output_Manager(
            mesh, nothing, atmo, 0, 600, [:not_a_variable];
            filename = joinpath(dir, "invalid.nc"),
            model_mode = :Barotropic,
        )
        @test_throws ArgumentError Output_Manager(
            mesh, nothing, atmo, 0, 600, [:u];
            filename = joinpath(dir, "invalid.nc"),
            model_mode = :Barotropic,
            output_interval = 0,
        )
    end
end

@testset "Pressure-level output" begin
    mesh = Spectral_Spherical_Mesh(1, 2, 8, 4, 2, 6.371e6)
    vert = Vert_Coordinate(
        8, 4, 2, "even_sigma", "simmons_and_burridge", "second_centered_wts",
    )
    atmo = Atmo_Data(
        "plev", 8, 4, 2, true, true, false, false, mesh.sinθ;
        radius = mesh.radius,
        grav = 10.0,
        rdgas = 300.0,
    )
    dyn = Dyn_Data("plev", 1, 2, 8, 4, 2)
    dyn.grid_t_c[:, :, 1] .= 250.0
    dyn.grid_t_c[:, :, 2] .= 290.0
    dyn.grid_u_c[:, :, 1] .= 10.0
    dyn.grid_u_c[:, :, 2] .= 20.0
    dyn.grid_z_full[:, :, 1] .= 8_000.0
    dyn.grid_z_full[:, :, 2] .= 1_000.0
    dyn.grid_ps_c .= 80_000.0
    requested = [:t, :u, :z]

    mktempdir() do dir
        filename = joinpath(dir, "primitive.nc")
        manager = Output_Manager(
            mesh, vert, atmo, 0, 600, requested;
            filename,
            model_mode = :PrimitiveEquation,
            do_plev_output = true,
            pressure_levels = [25_000.0, 75_000.0],
            output_interval = 3600,
        )
        @test requested == [:t, :u, :z]
        @test manager.active_symbols == [:t, :u, :z, :ps]
        Update_Output!(manager, dyn, 600)
        dyn.grid_t_c[:, :, 1] .= 260.0
        dyn.grid_t_c[:, :, 2] .= 300.0
        dyn.grid_u_c[:, :, 1] .= 12.0
        dyn.grid_u_c[:, :, 2] .= 22.0
        dyn.grid_ps_c .= 120_000.0
        Update_Output!(manager, dyn, 1200)
        Finalize_Output!(manager)

        pfull_top = 50_000.0 / exp(1.0)
        pfull_bottom = exp(
            (100_000.0 * log(100_000.0) - 50_000.0 * log(50_000.0)) /
            50_000.0 - 1.0,
        )
        weight_25k = (log(25_000.0) - log(pfull_bottom)) /
                     (log(pfull_top) - log(pfull_bottom))
        expected_25k = 295.0 + weight_25k * (255.0 - 295.0)

        NCDataset(joinpath(dir, "primitive_t0_plev.nc"), "r") do ds
            @test ds["plev"].var[:] == [25_000.0, 75_000.0]
            temperature = ds["ta"].var[:, :, :, 1]
            @test all(temperature[:, :, 1] .≈ expected_25k)
            @test all(isnan, temperature[:, :, 2])
        end

        # Isca-style native files retain approximate reference axes together
        # with the interface coefficients needed for exact reconstruction.
        NCDataset(joinpath(dir, "primitive_t0.nc"), "r") do ds
            @test ds["pfull"].var[:] ≈ [pfull_top, pfull_bottom] .* 0.01
            @test ds["phalf"].var[:] == [0.0, 500.0, 1000.0]
            @test ds["pk"].var[:] == vert.ak
            @test ds["bk"].var[:] == vert.bk
            @test ds.attrib["gravity"] == atmo.grav
            @test ds.attrib["dry_air_gas_constant"] == atmo.rdgas
            @test !haskey(ds, "hyam")
            @test !haskey(ds, "hybm")
        end

        offline_path = joinpath(dir, "primitive_offline_plev.nc")
        Interpolate_File(
            joinpath(dir, "primitive_t0.nc"), offline_path,
            [25_000.0, 75_000.0]; var_names=[:t, :u, :z],
        )
        NCDataset(offline_path, "r") do ds
            @test haskey(ds, "ta")
            @test haskey(ds, "ua")
            @test haskey(ds, "zg")
            @test !haskey(ds, "t")
            @test !haskey(ds, "u")
            temperature = ds["ta"].var[:, :, :, 1]
            @test all(temperature[:, :, 1] .≈ expected_25k)
            @test all(isnan, temperature[:, :, 2])
            wind = ds["ua"].var[:, :, :, 1]
            expected_wind = 21.0 + weight_25k * (11.0 - 21.0)
            @test all(wind[:, :, 1] .≈ expected_wind)
            @test all(isnan, wind[:, :, 2])

            NCDataset(joinpath(dir, "primitive_t0_plev.nc"), "r") do online_ds
                @test isequal(ds["zg"].var[:], online_ds["zg"].var[:])
            end
        end

        # Preserve Isca's current ordering: average the native record first,
        # then interpolate that record to pressure levels.
        function interpolate_sample(t_top, t_bottom, ps)
            sample_top = 0.5 * ps / exp(1.0)
            sample_bottom = ps * pfull_bottom / 100_000.0
            weight = (log(25_000.0) - log(sample_bottom)) /
                     (log(sample_top) - log(sample_bottom))
            return t_bottom + weight * (t_top - t_bottom)
        end
        interpolate_then_average = 0.5 * (
            interpolate_sample(250.0, 290.0, 80_000.0) +
            interpolate_sample(260.0, 300.0, 120_000.0)
        )
        @test !isapprox(expected_25k, interpolate_then_average; atol=1.0e-8)
    end

    @test_throws DimensionMismatch Vert_Coordinate(
        8, 4, 10, "simmons_and_burridge", "simmons_and_burridge",
        "second_centered_wts",
    )
end

@testset "Isca-compatible full-level pressure reconstruction" begin
    pressure = zeros(Float64, 1, 1, 2)
    ps = fill(100_000.0, 1, 1)

    JGCM.Vertical_Interpolation_Module.Compute_Pressure_Grid!(
        pressure, [0.0, 0.0, 0.0], [0.0, 0.5, 1.0], ps,
    )
    @test pressure[1, 1, 1] ≈ 50_000.0 / exp(1.0)
    @test pressure[1, 1, 2] ≈ exp(
        (100_000.0 * log(100_000.0) - 50_000.0 * log(50_000.0)) /
        50_000.0 - 1.0,
    )

    # Finite-top behavior is inferred from the coefficients themselves.
    JGCM.Vertical_Interpolation_Module.Compute_Pressure_Grid!(
        pressure, [200.0, 50_000.0, 100_000.0], zeros(3), ps,
    )
    expected_finite_top = exp(
        (50_000.0 * log(50_000.0) - 200.0 * log(200.0)) / 49_800.0 - 1.0,
    )
    @test pressure[1, 1, 1] ≈ expected_finite_top
    @test pressure[1, 1, 1] != 50_000.0 / exp(1.0)

    # The dynamics and output reconstruction must infer the same finite top
    # directly from the prescribed Simmons--Burridge interface coefficients.
    finite_top_vert = Vert_Coordinate(
        1, 1, 20, "simmons_and_burridge", "simmons_and_burridge",
        "second_centered_wts",
    )
    @test !finite_top_vert.zero_top
    dynamics_p_half = zeros(Float64, 1, 1, 21)
    dynamics_delta_p = zeros(Float64, 1, 1, 20)
    dynamics_log_p_half = zeros(Float64, 1, 1, 21)
    dynamics_p_full = zeros(Float64, 1, 1, 20)
    dynamics_log_p_full = zeros(Float64, 1, 1, 20)
    JGCM.Press_And_Geopot_Module.Pressure_Variables!(
        finite_top_vert,
        reshape(ps, 1, 1, 1),
        dynamics_p_half,
        dynamics_delta_p,
        dynamics_log_p_half,
        dynamics_p_full,
        dynamics_log_p_full,
    )
    output_p_full = similar(dynamics_p_full)
    JGCM.Vertical_Interpolation_Module.Compute_Pressure_Grid!(
        output_p_full, finite_top_vert.ak, finite_top_vert.bk, ps,
    )
    @test dynamics_p_full ≈ output_p_full
    @test dynamics_log_p_half[1, 1, 1] ≈ log(finite_top_vert.ak[1])

    finite_top_atmo = Atmo_Data(
        "finite top", 1, 1, 20, false, false, false, false, [0.0];
        radius=6.371e6,
    )
    virtual_temperature = fill(250.0, 1, 1, 20)
    geopotential_full = zeros(Float64, 1, 1, 20)
    geopotential_half = zeros(Float64, 1, 1, 21)
    JGCM.Press_And_Geopot_Module.Compute_Geopotential!(
        finite_top_vert,
        finite_top_atmo,
        dynamics_log_p_half,
        dynamics_log_p_full,
        virtual_temperature,
        zeros(Float64, 1, 1, 1),
        geopotential_full,
        geopotential_half,
    )
    for k in 1:20
        @test geopotential_half[1, 1, k] - geopotential_half[1, 1, k+1] ≈
              finite_top_atmo.rdgas * virtual_temperature[1, 1, k] *
              (dynamics_log_p_half[1, 1, k+1] - dynamics_log_p_half[1, 1, k])
    end

    # :z consistently exposes geometric height; geopotential remains an
    # internal dynamical field with different units.
    height_dyn = Dyn_Data("height mapping", 1, 2, 1, 1, 2)
    height_dyn.grid_z_full .= 1_500.0
    height_dyn.grid_geopot_full .= 14_700.0
    primitive_map = JGCM.Variable_Mappings_Module.Get_Dyn_Var_Map(
        height_dyn, Val(:PrimitiveEquation),
    )
    @test primitive_map[:z] === height_dyn.grid_z_full
    @test all(primitive_map[:z] .== 1_500.0)

    # Geopotential height follows Isca's hypsometric interpolation and remains
    # in meters (R/g), rather than accidentally producing geopotential units.
    atmo = Atmo_Data(
        "height interpolation", 1, 1, 2, false, false, false, false, [0.0];
        radius=6.371e6,
    )
    height = reshape([1000.0 + log(4.0) * 280.0 * atmo.rdgas / atmo.grav, 1000.0], 1, 1, 2)
    temperature = fill(280.0, 1, 1, 2)
    pfull = reshape([20_000.0, 80_000.0], 1, 1, 2)
    interpolated = zeros(Float64, 1, 1, 2)
    JGCM.Vertical_Interpolation_Module.Interpolate_Field!(
        interpolated, height, pfull, ps, log.([40_000.0, 90_000.0]),
        :z, atmo, temperature,
    )
    @test interpolated[1, 1, 1] ≈ 1000.0 + log(2.0) * 280.0 * atmo.rdgas / atmo.grav
    @test isnan(interpolated[1, 1, 2])
end

@testset "One-step driver checkpoint" begin
    mktempdir() do dir
        output_path = joinpath(dir, "nested", "run")
        config = Model_Config(
            name = "one_step",
            model_type = :Barotropic,
            num_fourier = 5,
            nθ = 12,
            nd = 1,
            radius = 6.371e6,
            omega = 7.292e-5,
            grav = 9.8,
            vert_coord_option = "even_sigma",
            vert_difference_option = "none",
            vert_ref_level_option = "none",
            Δt = 600,
            end_time = 600,
            day_to_sec = 86_400,
            damping_order = 4,
            damping_coef = 1.0e-4,
            robert_coef = 0.04,
            implicit_coef = 0.0,
            saving_frequency = 600,
            moisture_processes = false,
            initial_condition = :Barotropic_Jet,
            output_path = output_path,
            output_filename = joinpath(output_path, "output.nc"),
            logger = joinpath(output_path, "run.log"),
            vars_to_output = [:u],
            output_interval = 600,
            physics_params = Dict{String,Any}(),
        )

        JGCM_Simulate(config)
        @test isfile(joinpath(output_path, "output_t0.nc"))
        @test isfile(joinpath(output_path, "restart", "restart_t600.jld2"))
    end
end
