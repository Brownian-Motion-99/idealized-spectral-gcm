using FFTW
using LinearAlgebra
using Random
using Test
using JGCM

const Semi_Implicit_Kernels = JGCM.Semi_Implicit_Module

function direct_grid_to_spherical(mesh, grid, qwg)
    nf, ns = mesh.num_fourier, mesh.num_spherical
    fourier = fft(complex.(grid), (1,)) / mesh.nλ
    spectral = zeros(ComplexF64, nf + 1, ns + 1, size(grid, 3))
    for k = axes(grid, 3), m = 0:nf, n = m:ns-1
        spectral[m+1, n+1, k] =
            0.5 * sum(
                @view(fourier[m+1, :, k]) .*
                @view(qwg[m+1, n+1, :]),
            )
    end
    return spectral
end

function direct_spherical_to_grid(mesh, spectral, qnm)
    nf, ns = mesh.num_fourier, mesh.num_spherical
    fourier = zeros(ComplexF64, mesh.nλ, mesh.nθ, size(spectral, 3))
    for k = axes(spectral, 3), j = 1:mesh.nθ, m = 0:nf
        for n = m:ns
            fourier[m+1, j, k] +=
                spectral[m+1, n+1, k] * qnm[m+1, n+1, j]
        end
    end
    @views fourier[1, :, :] ./= 2
    return 2mesh.nλ .* real.(ifft(fourier, (1,)))
end

function allocating_mode_solve!(spectral, wave_matrix, num_wavenumbers)
    nf, ns, _ = size(spectral)
    for m = 0:nf-1, n = m:ns-1
        if n <= num_wavenumbers
            spectral[m+1, n+1, :] .=
                wave_matrix[:, :, n+1] * spectral[m+1, n+1, :]
        end
    end
    return nothing
end

function allocation_free_mode_solve!(
    spectral,
    wave_matrix,
    num_wavenumbers,
    work,
)
    nf, ns, _ = size(spectral)
    for m = 0:nf-1, n = m:ns-1
        if n <= num_wavenumbers
            mode = @view spectral[m+1, n+1, :]
            copyto!(work, mode)
            mul!(mode, @view(wave_matrix[:, :, n+1]), work)
        end
    end
    return nothing
end

@testset "Independent spectral transform validation" begin
    rng = MersenneTwister(2026)
    mesh = Spectral_Spherical_Mesh(3, 4, 24, 12, 2, 1.0)
    qnm, _ = Compute_Legendre(
        mesh.num_fourier,
        mesh.num_spherical,
        mesh.sinθ,
        mesh.nθ,
    )
    qwg = qnm .* reshape(mesh.wts, 1, 1, :)

    grid = randn(rng, mesh.nλ, mesh.nθ, mesh.nd)
    spectral = zeros(ComplexF64, 4, 5, mesh.nd)
    Trans_Grid_To_Spherical!(mesh, grid, spectral)
    @test spectral ≈
          direct_grid_to_spherical(mesh, grid, qwg) rtol = 2e-14 atol = 2e-14

    spectral_input = zeros(ComplexF64, 4, 5, mesh.nd)
    for k = 1:mesh.nd, m = 0:mesh.num_fourier, n = m:mesh.num_spherical
        spectral_input[m+1, n+1, k] = randn(rng) + im * randn(rng)
    end
    transformed_grid = similar(grid)
    Trans_Spherical_To_Grid!(mesh, spectral_input, transformed_grid)
    @test transformed_grid ≈
          direct_spherical_to_grid(mesh, spectral_input, qnm) rtol = 2e-14 atol = 2e-14
end

@testset "Batched Helmholtz solve equivalence and allocation" begin
    rng = MersenneTwister(2027)
    nf, ns, nd = 8, 9, 5
    num_wavenumbers = ns - 2
    wave_matrix = randn(rng, ComplexF64, nd, nd, num_wavenumbers + 1)
    initial = randn(rng, ComplexF64, nf, ns, nd)

    allocating = copy(initial)
    per_mode = copy(initial)
    batched = copy(initial)
    mode_work = zeros(ComplexF64, nd)
    batch_work = zeros(ComplexF64, nf, ns, nd)

    allocating_mode_solve!(allocating, wave_matrix, num_wavenumbers)
    allocation_free_mode_solve!(per_mode, wave_matrix, num_wavenumbers, mode_work)
    Semi_Implicit_Kernels.Helmholtz_Solve!(
        batched,
        wave_matrix,
        num_wavenumbers,
        batch_work,
    )

    @test batched ≈ allocating rtol = 5e-14 atol = 5e-14
    @test batched ≈ per_mode rtol = 5e-14 atol = 5e-14

    copyto!(batched, initial)
    Semi_Implicit_Kernels.Helmholtz_Solve!(
        batched,
        wave_matrix,
        num_wavenumbers,
        batch_work,
    )
    copyto!(batched, initial)
    allocations = @allocated Semi_Implicit_Kernels.Helmholtz_Solve!(
        batched,
        wave_matrix,
        num_wavenumbers,
        batch_work,
    )
    # BLAS may perform a few hundred bytes of thread-local bookkeeping on
    # some Julia/library combinations. Guard against array-sized regressions
    # without making the test runtime-version dependent.
    @test allocations <= 1024
end
