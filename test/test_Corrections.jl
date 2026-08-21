@testset "Water correction keeps grid and spectral states consistent" begin
    num_fourier, num_spherical = 3, 4
    nλ, nθ, nd = 64, 32, 2

    mesh = Spectral_Spherical_Mesh(
        num_fourier, num_spherical, nλ, nθ, nd, 6.371e6,
    )
    vert_coord = Vert_Coordinate(
        nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts",
    )
    atmo = Atmo_Data(
        "water correction test", nλ, nθ, nd,
        false, false, true, true,
        mesh.sinθ;
        radius=6.371e6,
    )
    data = Dyn_Data(
        "water correction test", num_fourier, num_spherical, nλ, nθ, nd, 1,
    )

    data.grid_ps_n .= 100_000.0
    for k in 1:nd, j in 1:nθ, i in 1:nλ
        data.grid_q_n[i, j, k] = 0.004 + 0.002 * cos(2pi * (i - 1) / nλ) *
                                mesh.sinθ[j] + 0.001 * (k - 1)
    end
    data.grid_q_n[1, 1, 1] = -0.003

    target_water = 0.8 * JGCM.Vert_Coordinate_Module.Mass_Weighted_Global_Integral(
        vert_coord, mesh, atmo, max.(data.grid_q_n, 0.0), data.grid_ps_n,
    )

    JGCM.Spectral_Dynamics_Module.Compute_Corrections!(
        mesh, vert_coord, atmo,
        0.0, data.grid_ps_n, data.spe_lnps_n,
        0.0, data.grid_energy_full,
        data.grid_u_n, data.grid_v_n, data.grid_t_n, data.spe_t_n,
        target_water, data.grid_q_n, data.spe_q_n,
    )

    reconstructed = similar(data.grid_q_n)
    Trans_Spherical_To_Grid!(mesh, data.spe_q_n, reconstructed)
    corrected_water = JGCM.Vert_Coordinate_Module.Mass_Weighted_Global_Integral(
        vert_coord, mesh, atmo, data.grid_q_n, data.grid_ps_n,
    )

    @test data.grid_q_n ≈ reconstructed atol=1e-14 rtol=1e-13
    @test corrected_water ≈ target_water rtol=1e-13

    JGCM.Dyn_Data_Module.Time_Advance!(data)
    Trans_Spherical_To_Grid!(mesh, data.spe_q_c, reconstructed)
    @test data.grid_q_c ≈ reconstructed atol=1e-14 rtol=1e-13
end

@testset "Zero-water correction is finite" begin
    num_fourier, num_spherical = 3, 4
    nλ, nθ, nd = 64, 32, 2

    mesh = Spectral_Spherical_Mesh(
        num_fourier, num_spherical, nλ, nθ, nd, 6.371e6,
    )
    vert_coord = Vert_Coordinate(
        nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts",
    )
    atmo = Atmo_Data(
        "zero-water correction test", nλ, nθ, nd,
        false, false, true, true,
        mesh.sinθ;
        radius=6.371e6,
    )
    data = Dyn_Data(
        "zero-water correction test", num_fourier, num_spherical, nλ, nθ, nd, 1,
    )
    data.grid_ps_n .= 100_000.0
    data.grid_q_n .= 0.001
    data.spe_q_n .= 1.0 + 1.0im

    JGCM.Spectral_Dynamics_Module.Compute_Corrections!(
        mesh, vert_coord, atmo,
        0.0, data.grid_ps_n, data.spe_lnps_n,
        0.0, data.grid_energy_full,
        data.grid_u_n, data.grid_v_n, data.grid_t_n, data.spe_t_n,
        0.0, data.grid_q_n, data.spe_q_n,
    )

    @test all(iszero, data.grid_q_n)
    @test all(iszero, data.spe_q_n)
    @test all(isfinite, data.grid_q_n)
    @test all(isfinite, data.spe_q_n)
end
