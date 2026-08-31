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

    # Isca-style large-Courant zonal transport crosses whole cells without
    # subcycling when the flow is coherent (neighboring departure faces do
    # not cross). The transverse half-step uses the same departure-point
    # interpolation for either wind direction.
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
    @test steps.horizontal_substeps == 1
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
    # and deformation-only safety subcycling.
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

@testset "Tracer roundoff undershoots" begin
    cleanup! = JGCM.Grid_Tracer_Transport_Module._remove_roundoff_undershoots!
    tracer = [0.01, -4.0e-14]
    cleanup!(tracer, "test")
    @test tracer == [0.01, 0.0]
    @test_throws ErrorException cleanup!([0.01, -1.0e-8], "test")
end
