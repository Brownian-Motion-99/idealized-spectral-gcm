using Test
using JGCM

function column_pressure_budgets(vert, ps, q)
    nλ, nθ, nd = size(q)
    dry = zeros(nλ, nθ)
    water = zeros(nλ, nθ)
    layer_water = similar(q)
    for k in 1:nd
        @views Δp = @. vert.Δak[k] + vert.Δbk[k] * ps[:, :, 1]
        @views @. layer_water[:, :, k] = q[:, :, k] * Δp
        @views @. water += layer_water[:, :, k]
        @views @. dry += (1.0 - q[:, :, k]) * Δp
    end
    return dry, water, layer_water
end

@testset "Dry-air pressure adjustment" begin
    for coordinate in ("even_sigma", "hybrid")
        nλ, nθ, nd = 2, 2, 6
        vert = Vert_Coordinate(
            nλ,
            nθ,
            nd,
            coordinate,
            "simmons_and_burridge",
            "second_centered_wts",
        )
        ps = reshape([100_000.0, 98_000.0, 96_000.0, 94_000.0], nλ, nθ, 1)
        ps_before = copy(ps)
        q_before = zeros(nλ, nθ, nd)
        for k in 1:nd
            @views q_before[:, :, k] .= 0.0015 * k
        end
        q_physics = copy(q_before)
        q_physics[:, :, 2] .-= 4.0e-4
        q_physics[:, :, 5] .+= 1.5e-4
        q_adjusted = copy(q_physics)

        Δp_before = similar(q_before)
        for k in 1:nd
            @views @. Δp_before[:, :, k] =
                vert.Δak[k] + vert.Δbk[k] * ps_before[:, :, 1]
        end
        dry_before, _, _ = column_pressure_budgets(vert, ps_before, q_before)
        _, water_physics, layer_water_physics =
            column_pressure_budgets(vert, ps_before, q_physics)

        Dry_Air_Adjustment!(vert, ps, q_before, q_adjusted, Δp_before)

        dry_after, water_after, layer_water_after =
            column_pressure_budgets(vert, ps, q_adjusted)
        expected_pressure_change = dropdims(
            sum((q_physics .- q_before) .* Δp_before; dims = 3);
            dims = 3,
        ) ./ sum(vert.Δbk)
        @test dropdims(ps .- ps_before; dims = 3) ≈ expected_pressure_change
        @test dry_after ≈ dry_before rtol = 2.0e-14 atol = 1.0e-9
        @test water_after ≈ water_physics rtol = 2.0e-14 atol = 1.0e-10
        @test layer_water_after ≈ layer_water_physics rtol = 2.0e-14 atol = 1.0e-10
        @test all(0.0 .<= q_adjusted .< 1.0)
    end

    nλ, nθ, nd = 1, 1, 4
    vert = Vert_Coordinate(
        nλ,
        nθ,
        nd,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )
    ps = fill(100_000.0, nλ, nθ, 1)
    q = reshape([0.001, 0.003, 0.006, 0.012], nλ, nθ, nd)
    Δp = fill(25_000.0, nλ, nθ, nd)
    ps_identity, q_identity = copy(ps), copy(q)
    Dry_Air_Adjustment!(vert, ps_identity, q, q_identity, Δp)
    @test ps_identity == ps
    @test q_identity == q

    invalid_q = copy(q)
    invalid_q[1] = 1.0
    @test_throws DomainError Dry_Air_Adjustment!(vert, ps, q, invalid_q, Δp)
    invalid_Δp = copy(Δp)
    invalid_Δp[1] = 0.0
    @test_throws DomainError Dry_Air_Adjustment!(vert, ps, q, copy(q), invalid_Δp)
end
