using Test
using JGCM
using LinearAlgebra

@testset "Gaussian quadrature and Legendre basis" begin
    num_fourier, nθ = 7, 16
    num_spherical = num_fourier + 1
    sinθ, wts = Compute_Gaussian(nθ)
    qnm, dqnm = Compute_Legendre(num_fourier, num_spherical, sinθ, nθ)

    q44 = sqrt(35) / 4 * (1 .- sinθ .^ 2) .^ (3 / 2)
    q34 = sqrt(105 / 8) * (sinθ .- sinθ .^ 3)
    q13 = sqrt(5) / 2 * (3sinθ .^ 2 .- 1)
    q14 = sqrt(7) / 2 * (5sinθ .^ 3 .- 3sinθ)
    q24 = sqrt(21) / 4 * (5sinθ .^ 2 .- 1) .* sqrt.(1 .- sinθ .^ 2)


    dq44 = sqrt(35) / 4 * (1 .- sinθ .^ 2) .^ (1 / 2) * (3 / 2) .* (-2sinθ)
    dq34 = sqrt(105 / 8) * (1 .- 3 * sinθ .^ 2)
    dq13 = sqrt(5) / 2 * (6sinθ)
    dq14 = sqrt(7) / 2 * (15sinθ .^ 2 .- 3)
    dq24 =
        sqrt(21) / 4 * (10 * sinθ) .* sqrt.(1 .- sinθ .^ 2) +
        sqrt(21) / 4 * (5sinθ .^ 2 .- 1) ./ sqrt.(1 .- sinθ .^ 2) * 0.5 .* (-2sinθ)

    for (actual, expected) in (
        (qnm[4, 4, :], q44), (qnm[3, 4, :], q34),
        (qnm[1, 3, :], q13), (qnm[1, 4, :], q14),
        (qnm[2, 4, :], q24), (dqnm[4, 4, :], dq44),
        (dqnm[3, 4, :], dq34), (dqnm[1, 3, :], dq13),
        (dqnm[1, 4, :], dq14), (dqnm[2, 4, :], dq24),
    )
        @test actual ≈ expected atol = 1.0e-12
    end

    #check 1/2∫_{-1}^{1}P_{n,m} P_{l,m} = δ_{n,l}

    max_orthogonality_error = 0.0
    for m = 1:num_fourier+1
        for l = m:num_spherical+1
            for n = m:num_spherical+1
                projection = 0.5 * sum(qnm[m, n, :] .* qnm[m, l, :] .* wts)
                max_orthogonality_error = max(
                    max_orthogonality_error,
                    abs(projection - Float64(n == l)),
                )
            end
        end
    end
    @test max_orthogonality_error < 1.0e-10
end
