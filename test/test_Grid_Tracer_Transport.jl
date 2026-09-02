using Test
using JGCM

@testset "Grid tracer transport" begin
    nλ, nθ, nd = 32, 16, 5
    radius = 6.371e6
    mesh = Spectral_Spherical_Mesh(7, 8, nλ, nθ, nd, radius)
    workspace = Grid_Tracer_Workspace(nλ, nθ, nd)
    u = zeros(nλ, nθ, nd)
    v = zeros(nλ, nθ, nd)
    Δp = fill(20_000.0, nλ, nθ, nd)
    M = zeros(nλ, nθ, nd + 1)
    q0 = fill(0.01, nλ, nθ, nd)
    q1 = similar(q0)

    # Both advective-form operators must cancel exactly for a constant tracer,
    # including divergent horizontal flow and nonzero vertical mass flux.
    for k = 1:nd, j = 1:nθ, i = 1:nλ
        u[i, j, k] = 25sin(mesh.λc[i]) * mesh.cosθ[j]
        v[i, j, k] = 8cos(mesh.λc[i]) * mesh.sinθ[j]
    end
    M[:, :, 2:nd] .= 12.0
    steps = Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, 600.0)
    @test q1 ≈ q0 atol = 2e-17 rtol = 2e-15
    @test steps.horizontal_substeps >= 1
    @test steps.vertical_substeps >= 1
    @test @allocated(Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, 600.0)) == 0

    # Zonal transport is periodic, monotone, and conservative for this
    # nondivergent flow.
    fill!(v, 0.0)
    fill!(M, 0.0)
    u .= 40.0
    for k = 1:nd, j = 1:nθ, i = 1:nλ
        q0[i, j, k] = 0.002 + 0.01exp(-8sin(mesh.λc[i] - π)^2) * mesh.cosθ[j]^2
    end
    initial_mean = sum(q0 .* reshape(mesh.wts, 1, nθ, 1))
    Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, 1200.0)
    final_mean = sum(q1 .* reshape(mesh.wts, 1, nθ, 1))
    @test minimum(q1) >= 0.0
    @test minimum(q1) >= minimum(q0) - 2e-15
    @test maximum(q1) <= maximum(q0) + 2e-15
    @test final_mean ≈ initial_mean rtol = 2e-14

    # The swept reconstruction itself crosses whole cells, and the transverse
    # half-step uses the same departure-point interpolation for either wind
    # direction. The complete advective-form update subcycles separately to
    # keep its donor-cell baseline positive.
    transverse! = JGCM.Grid_Tracer_Transport_Module._transverse_states!
    for k = 1:nd, j = 1:nθ, i = 1:nλ
        q0[i, j, k] = i
        u[i, j, k] = 2.5 * radius * mesh.cosθ[j] * (2π / nλ) / 600.0
    end
    transverse!(workspace, mesh, q0, u, v, 600.0)
    i = 10
    @test workspace.q_half_lambda[i, 8, 3] ≈ 0.25q0[i-2, 8, 3] + 0.75q0[i-1, 8, 3]
    u .*= -1
    transverse!(workspace, mesh, q0, u, v, 600.0)
    @test workspace.q_half_lambda[i, 8, 3] ≈ 0.75q0[i+1, 8, 3] + 0.25q0[i+2, 8, 3]

    for k = 1:nd, j = 1:nθ, i = 1:nλ
        q0[i, j, k] = 0.002 + 0.01exp(-8sin(mesh.λc[i] - π)^2) * mesh.cosθ[j]^2
        u[i, j, k] = 2.4 * radius * mesh.cosθ[j] * (2π / nλ) / 600.0
    end
    initial_mean = sum(q0 .* reshape(mesh.wts, 1, nθ, 1))
    steps = Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, 600.0)
    final_mean = sum(q1 .* reshape(mesh.wts, 1, nθ, 1))
    @test steps.horizontal_substeps > 1
    @test minimum(q1) >= 0.0
    @test final_mean ≈ initial_mean rtol = 3e-14

    # Multi-layer PPM weights fully crossed nonuniform layers by their mass,
    # then reconstructs only the terminal layer fraction.
    ppm! = JGCM.Grid_Tracer_Transport_Module._ppm_reconstruction!
    swept = JGCM.Grid_Tracer_Transport_Module._vertical_swept_averages
    Δp_nonuniform = similar(Δp)
    for k = 1:nd
        Δp_nonuniform[:, :, k] .= 10_000k
        q0[:, :, k] .= 0.001k
    end
    ppm!(workspace, q0, Δp_nonuniform)
    low_average, high_average = swept(workspace, q0, Δp_nonuniform, 1, 1, 4, 55_000.0)
    expected = (30_000q0[1, 1, 3] + 20_000q0[1, 1, 2] + 5_000q0[1, 1, 1]) / 55_000
    @test low_average ≈ expected
    @test high_average ≈ expected

    # Exercise both vertical donor directions, boundary-adjacent PPM cells,
    # and incoming-mass safety subcycling.
    fill!(u, 0.0)
    q0 .= 0.001
    q0[:, :, 3] .= 0.018
    M[:, :, 2] .= -35.0
    M[:, :, 3] .= 35.0
    M[:, :, 4] .= -35.0
    M[:, :, 5] .= 35.0
    steps = Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, 600.0)
    @test steps.vertical_substeps == 3
    @test all(isfinite, q1)
    @test minimum(q1) >= 0.0
    @test maximum(q1) <= maximum(q0) + 2e-15
end


@testset "High-Courant dry-front positivity" begin
    # A gradual acceleration followed by sharp zonal convergence. The
    # deformation-only estimator selected one substep here; even a two-sided
    # deformation estimator selected only two and left the donor baseline
    # negative because the downstream face still swept the complete wet cell.
    nλ, nθ, nd = 32, 16, 1
    radius = 6.371e6
    dt = 600.0
    mesh = Spectral_Spherical_Mesh(7, 8, nλ, nθ, nd, radius)
    workspace = Grid_Tracer_Workspace(nλ, nθ, nd)
    u = zeros(nλ, nθ, nd)
    v = zeros(nλ, nθ, nd)
    Δp = fill(20_000.0, nλ, nθ, nd)
    M = zeros(nλ, nθ, nd + 1)
    Δλ = 2π / nλ
    for j = 1:nθ, i = 1:nλ
        center_courant = 4.0 + 2.0 * (i - 1) / (nλ - 1)
        u[i, j, 1] = center_courant * radius * mesh.cosθ[j] * Δλ / dt
    end

    q0 = zeros(nλ, nθ, nd)
    q0[1, 8, 1] = 1.0e-12
    q1 = similar(q0)
    steps = Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, dt)
    @test steps.horizontal_substeps > 2
    @test all(isfinite, q1)
    @test minimum(q1) >= 0.0
    @test maximum(q1) <= maximum(q0) + 2e-15

    # Constant tracers must remain constant under the same strongly varying
    # velocity field despite the additional positivity substeps.
    fill!(q0, 0.01)
    Advance_Grid_Tracer!(workspace, mesh, q1, q0, u, v, Δp, M, dt)
    @test q1 ≈ q0 atol = 2e-17 rtol = 3e-15

    # Vertical analogue: large coherent positive mass flux with a small local
    # convergence. A deformation-only bound leaves a multi-layer downstream
    # sweep and can empty the wet layer before applying the local correction.
    nλv, nθv, ndv = 4, 2, 20
    meshv = Spectral_Spherical_Mesh(1, 2, nλv, nθv, ndv, radius)
    workspacev = Grid_Tracer_Workspace(nλv, nθv, ndv)
    uv = zeros(nλv, nθv, ndv)
    vv = similar(uv)
    Δpv = ones(nλv, nθv, ndv)
    Mv = zeros(nλv, nθv, ndv + 1)
    for h = 2:10
        Mv[:, :, h] .= 5.0 * (h - 1) / 9
    end
    Mv[:, :, 11] .= 4.0
    for h = 12:20
        Mv[:, :, h] .= 4.0 * (21 - h) / 10
    end
    qv0 = zeros(nλv, nθv, ndv)
    qv0[:, :, 10] .= 1.0e-12
    qv1 = similar(qv0)
    steps = Advance_Grid_Tracer!(workspacev, meshv, qv1, qv0, uv, vv, Δpv, Mv, 1.0)
    @test steps.vertical_substeps > 2
    @test all(isfinite, qv1)
    @test minimum(qv1) >= 0.0

    fill!(qv0, 0.01)
    Advance_Grid_Tracer!(workspacev, meshv, qv1, qv0, uv, vv, Δpv, Mv, 1.0)
    @test qv1 ≈ qv0 atol = 2e-17 rtol = 3e-15
end

@testset "Tracer roundoff undershoots" begin
    cleanup! = JGCM.Grid_Tracer_Transport_Module._remove_roundoff_undershoots!
    tracer = [0.01, -4.0e-14]
    cleanup!(tracer, "test")
    @test tracer == [0.01, 0.0]
    @test_throws ErrorException cleanup!([0.01, -1.0e-8], "test")
end
