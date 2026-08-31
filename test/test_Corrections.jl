@testset "Grid-only water correction restores its requested integral" begin
    num_fourier, num_spherical = 3, 4
    nλ, nθ, nd = 64, 32, 2

    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nλ, nθ, nd, 6.371e6)
    vert_coord = Vert_Coordinate(
        nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts",
    )
    atmo = Atmo_Data(
        "water correction test", nλ, nθ, nd,
        false, false, true, true,
        mesh.sinθ;
        radius = 6.371e6,
    )
    data = Dyn_Data(
        "water correction test", num_fourier, num_spherical, nλ, nθ, nd, 1,
    )

    data.grid_ps_n .= 100_000.0
    for k = 1:nd, j = 1:nθ, i = 1:nλ
        data.grid_q_n[i, j, k] =
            0.004 + 0.002 * cos(2pi * (i - 1) / nλ) * mesh.sinθ[j] + 0.001 * (k - 1)
    end
    target_water = 0.8 * JGCM.Vert_Coordinate_Module.Mass_Weighted_Global_Integral(
        vert_coord, mesh, atmo, data.grid_q_n, data.grid_ps_n,
    )

    JGCM.Spectral_Dynamics_Module.Compute_Corrections!(
        mesh, vert_coord, atmo,
        0.0, data.grid_ps_n, data.spe_lnps_n,
        0.0, data.grid_energy_full,
        data.grid_u_n, data.grid_v_n, data.grid_t_n, data.spe_t_n,
        target_water, data.grid_q_n,
    )

    corrected_water = JGCM.Vert_Coordinate_Module.Mass_Weighted_Global_Integral(
        vert_coord, mesh, atmo, data.grid_q_n, data.grid_ps_n,
    )

    @test minimum(data.grid_q_n) >= -1e-14
    @test corrected_water ≈ target_water rtol = 1e-13

    JGCM.Dyn_Data_Module.Time_Advance!(data)
    @test data.grid_q_c == data.grid_q_n
end

@testset "Filtered current humidity uses the supplied current pressure" begin
    nλ, nθ, nd = 16, 8, 3
    mesh = Spectral_Spherical_Mesh(3, 4, nλ, nθ, nd, 6.371e6)
    vert = Vert_Coordinate(
        nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts",
    )
    atmo = Atmo_Data(
        "filtered humidity", nλ, nθ, nd, false, false, true, true, mesh.sinθ;
        radius = 6.371e6,
    )
    integrator = JGCM.Time_Integrator_Module.Filtered_Leapfrog(
        0.04, 4, 1.0e-4, mesh.laplacian_eig, 0.5, 600, false, 0, 1200,
    )
    qp = fill(0.004, nλ, nθ, nd)
    qc = fill(0.008, nλ, nθ, nd)
    qn = fill(0.012, nλ, nθ, nd)
    ps_c = zeros(nλ, nθ, 1)
    for j = 1:nθ
        ps_c[:, j, 1] .= 90_000.0 + 8_000.0 * mesh.sinθ[j]
    end
    water_before = JGCM.Vert_Coordinate_Module.Mass_Weighted_Global_Integral(
        vert, mesh, atmo, qc, ps_c,
    )
    JGCM.Spectral_Dynamics_Module._filter_grid_humidity!(
        integrator, mesh, vert, atmo, qp, qc, qn, ps_c,
    )
    water_after = JGCM.Vert_Coordinate_Module.Mass_Weighted_Global_Integral(
        vert, mesh, atmo, qc, ps_c,
    )
    @test water_after ≈ water_before rtol = 2e-15
    @test minimum(qc) >= 0.0
end

@testset "Zero-water correction is finite" begin
    num_fourier, num_spherical = 3, 4
    nλ, nθ, nd = 64, 32, 2

    mesh = Spectral_Spherical_Mesh(num_fourier, num_spherical, nλ, nθ, nd, 6.371e6)
    vert_coord = Vert_Coordinate(
        nλ, nθ, nd, "even_sigma", "simmons_and_burridge", "second_centered_wts",
    )
    atmo = Atmo_Data(
        "zero-water correction test", nλ, nθ, nd,
        false, false, true, true,
        mesh.sinθ;
        radius = 6.371e6,
    )
    data = Dyn_Data(
        "zero-water correction test", num_fourier, num_spherical, nλ, nθ, nd, 1,
    )
    data.grid_ps_n .= 100_000.0
    data.grid_q_n .= 0.001

    JGCM.Spectral_Dynamics_Module.Compute_Corrections!(
        mesh, vert_coord, atmo,
        0.0, data.grid_ps_n, data.spe_lnps_n,
        0.0, data.grid_energy_full,
        data.grid_u_n, data.grid_v_n, data.grid_t_n, data.spe_t_n,
        0.0, data.grid_q_n,
    )

    @test all(iszero, data.grid_q_n)
    @test all(isfinite, data.grid_q_n)
end
