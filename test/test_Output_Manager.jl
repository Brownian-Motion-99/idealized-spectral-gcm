using Test
using NCDatasets
using JGCM

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
    )
    dyn = Dyn_Data("plev", 1, 2, 8, 4, 2)
    dyn.grid_t_c[:, :, 1] .= 250.0
    dyn.grid_t_c[:, :, 2] .= 290.0
    dyn.grid_ps_c .= 100_000.0
    requested = [:t]

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
        @test requested == [:t]
        @test manager.active_symbols == [:t, :ps]
        Update_Output!(manager, dyn, 600)
        dyn.grid_t_c[:, :, 1] .= 260.0
        dyn.grid_t_c[:, :, 2] .= 300.0
        Update_Output!(manager, dyn, 1200)
        Finalize_Output!(manager)

        NCDataset(joinpath(dir, "primitive_t0_plev.nc"), "r") do ds
            @test ds["plev"].var[:] == [25_000.0, 75_000.0]
            temperature = ds["ta"].var[:, :, :, 1]
            @test all(temperature[:, :, 1] .== 255.0)
            @test all(temperature[:, :, 2] .== 295.0)
        end
    end

    @test_throws DimensionMismatch Vert_Coordinate(
        8, 4, 10, "simmons_and_burridge", "simmons_and_burridge",
        "second_centered_wts",
    )
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
