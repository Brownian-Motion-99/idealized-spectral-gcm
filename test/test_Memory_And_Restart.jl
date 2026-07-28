using Test
using JLD2
using JGCM

mutable struct LegacyRestartState
    spe_vor_p::Array{ComplexF64,3}
    grid_u_c::Array{Float64,3}
end

@testset "Restart compatibility" begin
    mktempdir() do dir
        manager = Restart_Manager(dir, 600)
        source = Dyn_Data("source", 1, 2, 4, 2, 2, 0)
        source.spe_vor_p[1] = 3 + 4im
        source.grid_u_c[1] = 12.5
        Write_Restart_File(manager, source, 600)

        restored = Dyn_Data("restored", 1, 2, 4, 2, 2, 0)
        @test Load_Restart_File!(restored, joinpath(dir, "restart_t600.jld2")) == 600
        @test restored.spe_vor_p == source.spe_vor_p
        @test restored.grid_u_c == source.grid_u_c

        legacy_file = joinpath(dir, "legacy.jld2")
        legacy = LegacyRestartState(copy(source.spe_vor_p), copy(source.grid_u_c))
        jldsave(legacy_file; dyn_data_state = legacy, saved_time = 300)
        legacy_restored = Dyn_Data("legacy", 1, 2, 4, 2, 2, 0)
        @test Load_Restart_File!(legacy_restored, legacy_file) == 300
        @test legacy_restored.spe_vor_p == source.spe_vor_p
        @test legacy_restored.grid_u_c == source.grid_u_c
    end
end

@testset "Memory configuration" begin
    @test !hasfield(Dyn_Data, :spec_δu)
    @test !hasfield(Dyn_Data, :spec_δv)
    @test !hasfield(Dyn_Data, :grid_d_half1)
    @test !hasfield(Dyn_Data, :grid_d_half2)

    mesh = Spectral_Spherical_Mesh(3, 4, 64, 32, 2, 6.371e6)
    vert = Vert_Coordinate(
        64,
        32,
        2,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )
    atmo = Atmo_Data(
        "output",
        64,
        32,
        2,
        true,
        true,
        true,
        true,
        mesh.sinθ;
        radius = 6.371e6,
    )
    mktempdir() do dir
        output = Output_Manager(
            mesh,
            vert,
            atmo,
            0,
            600,
            [:t, :ps];
            filename = joinpath(dir, "output.nc"),
            do_plev_output = false,
            pressure_levels = [100_000.0, 50_000.0],
            output_interval = 600,
        )
        @test !isempty(output.acc_raw)
        Finalize_Output!(output)
    end
end

@testset "Reusable PBL workspace" begin
    mesh = Spectral_Spherical_Mesh(3, 4, 64, 32, 2, 6.371e6)
    atmo = Atmo_Data(
        "pbl",
        64,
        32,
        2,
        true,
        true,
        true,
        true,
        mesh.sinθ;
        radius = 6.371e6,
    )
    workspace = JGCM.Atmos_Param_Module.PBL_Workspace(64, 32, 2)
    p_half = fill(90_000.0, 64, 32, 3)
    p_full = fill(85_000.0, 64, 32, 2)
    ps = fill(100_000.0, 64, 32, 1)
    u = fill(3.0, 64, 32, 2)
    v = fill(4.0, 64, 32, 2)
    t = fill(280.0, 64, 32, 2)
    q = fill(0.01, 64, 32, 2)

    V_c, za, rho = JGCM.Atmos_Param_Module.Calculate_V_c_za_rho!(
        workspace,
        atmo,
        p_half,
        p_full,
        ps,
        u,
        v,
        t,
        q,
    )
    @test V_c === workspace.V_c
    @test za === workspace.za
    @test rho === workspace.rho
    @test all(V_c .== 5.0)
    @test all(isfinite, za)
    @test all(isfinite, rho)
end
