using JLD2

@testset "LRF state and kernel" begin
    nλ, nθ, nd = 2, 2, 2
    day_to_sec = 86400

    LW_q = zeros(Float64, nd, nd, nθ)
    LW_q[:, :, 1] .= [1.0 2.0; 3.0 4.0]
    LW_q[:, :, 2] .= [-1.0 0.5; 2.0 -2.0]

    ref_q = reshape(collect(1.0:(nλ*nθ*nd)), nλ, nθ, nd) .* 1.0e-3
    q = copy(ref_q)
    q[1, 1, :] .+= [1.0e-3, 2.0e-3]
    q[2, 1, :] .+= [-1.0e-3, 3.0e-3]
    q[1, 2, :] .+= [4.0e-3, -2.0e-3]

    state = LRF_State(LW_q, ref_q)
    tendency = fill(NaN, nλ, nθ, nd)
    LRF!(state, q, tendency, day_to_sec)

    @test tendency[1, 1, :] ≈ ([1.0 2.0; 3.0 4.0] * [1.0e-3, 2.0e-3]) ./ day_to_sec
    @test tendency[2, 1, :] ≈ ([1.0 2.0; 3.0 4.0] * [-1.0e-3, 3.0e-3]) ./ day_to_sec
    @test tendency[1, 2, :] ≈ ([-1.0 0.5; 2.0 -2.0] * [4.0e-3, -2.0e-3]) ./ day_to_sec
    @test tendency[2, 2, :] == zeros(nd)
    @test all(isfinite, tendency)

    fill!(tendency, 42.0)
    LRF!(state, ref_q, tendency, day_to_sec)
    @test tendency == zeros(size(tendency))

    @test_throws DimensionMismatch LRF_State(zeros(nd, nd + 1, nθ), ref_q)
    @test_throws DimensionMismatch LRF!(state, zeros(nλ + 1, nθ, nd), tendency, day_to_sec)
    @test_throws ArgumentError LRF!(state, ref_q, tendency, 0)

    mktempdir() do directory
        filepath = joinpath(directory, "lrf_test.jld2")
        JLD2.jldsave(filepath; LRF_LW_q = LW_q, ref_q = ref_q)

        loaded = Load_LRF_State(filepath, nλ, nθ, nd)
        @test loaded.LW_q == LW_q
        @test loaded.ref_q == ref_q

        @test_throws DimensionMismatch Load_LRF_State(filepath, nλ, nθ + 1, nd)
    end
end

@testset "Regularized logarithmic LRF" begin
    nλ, nθ, nd = 2, 3, 2
    day_to_sec = 86400
    q0 = 1.0e-8
    alpha = 0.8
    B = [2.0 -1.0; 0.5 3.0]
    taper = [1.0, 0.5, 0.0]
    reference_q = [0.0 2.0e-3 4.0e-3; 1.0e-4 5.0e-3 8.0e-3]
    chi_reference = log.(reference_q .+ q0)
    q = zeros(nλ, nθ, nd)
    for i = 1:nλ, j = 1:nθ
        q[i, j, :] .= reference_q[:, j]
    end

    state = Regularized_LRF_State(B, chi_reference, taper, q0, alpha, nλ)
    tendency = fill(NaN, size(q))
    LRF!(state, q, tendency, day_to_sec)
    @test tendency ≈ zeros(size(q)) atol = 2.0e-19
    @test all(isfinite, tendency)

    delta_chi = [0.1, 0.2]
    for j = 1:nθ
        q[1, j, :] .= exp.(chi_reference[:, j] .+ delta_chi) .- q0
    end
    LRF!(state, q, tendency, day_to_sec)
    expected_column = alpha .* (B * delta_chi) ./ day_to_sec
    @test tendency[1, 1, :] ≈ expected_column
    @test tendency[1, 2, :] ≈ 0.5 .* expected_column
    @test tendency[1, 3, :] == zeros(nd)

    invalid = copy(q)
    invalid[1, 1, 1] = -eps()
    @test_throws DomainError LRF!(state, invalid, tendency, day_to_sec)
    invalid[1, 1, 1] = NaN
    @test_throws DomainError LRF!(state, invalid, tendency, day_to_sec)
    @test_throws ArgumentError Regularized_LRF_State(B, chi_reference, [1.0, 1.1, 0.0], q0, alpha, nλ)
    @test_throws DimensionMismatch Regularized_LRF_State(B[:, 1:1], chi_reference, taper, q0, alpha, nλ)

    mktempdir() do directory
        filepath = joinpath(directory, "regularized_lrf_test.jld2")
        latitude = [-30.0, 0.0, 30.0]
        JLD2.jldsave(
            filepath;
            scheme = "regularized_log_tapered_v1",
            B_star = B,
            chi_reference = chi_reference,
            taper = taper,
            q0_kg_kg = q0,
            alpha = alpha,
            latitude = latitude,
        )
        loaded = Load_LRF_State(filepath, nλ, nθ, nd; latitude = latitude)
        @test loaded isa Regularized_LRF_State
        @test loaded.B == B
        @test loaded.chi_reference == chi_reference
        @test loaded.taper == taper
        @test loaded.q0 == q0
        @test loaded.alpha == alpha
        @test_throws ArgumentError Load_LRF_State(
            filepath, nλ, nθ, nd; latitude = reverse(latitude))
        @test_throws DimensionMismatch Load_LRF_State(filepath, nλ, nθ + 1, nd)
    end
end

@testset "Latitude-specific regularized logarithmic LRF" begin
    nλ, nθ, nd = 2, 3, 2
    day_to_sec = 86400
    q0 = 1.0e-8
    alpha = 0.7
    B_by_lat = zeros(nd, nd, nθ)
    B_by_lat[:, :, 1] .= [1.0 0.0; 0.0 2.0]
    B_by_lat[:, :, 2] .= [2.0 -1.0; 0.5 3.0]
    B_by_lat[:, :, 3] .= [-1.0 0.5; 2.0 -2.0]
    reference_q = [0.0 2.0e-3 4.0e-3; 1.0e-4 5.0e-3 8.0e-3]
    chi_reference = log.(reference_q .+ q0)
    q = zeros(nλ, nθ, nd)
    for i = 1:nλ, j = 1:nθ
        q[i, j, :] .= reference_q[:, j]
    end

    state = Latitude_LRF_State(B_by_lat, chi_reference, q0, alpha, nλ)
    tendency = fill(NaN, size(q))
    LRF!(state, q, tendency, day_to_sec)
    @test tendency ≈ zeros(size(q)) atol = 2.0e-19
    @test all(isfinite, tendency)

    delta_chi = [0.1, 0.2]
    for j = 1:nθ
        q[1, j, :] .= exp.(chi_reference[:, j] .+ delta_chi) .- q0
    end
    LRF!(state, q, tendency, day_to_sec)
    for j = 1:nθ
        @test tendency[1, j, :] ≈ alpha .* (B_by_lat[:, :, j] * delta_chi) ./ day_to_sec
    end
    @test tendency[2, :, :] ≈ zeros(nθ, nd) atol = 2.0e-19

    invalid = copy(q)
    invalid[1, 1, 1] = -eps()
    @test_throws DomainError LRF!(state, invalid, tendency, day_to_sec)
    @test_throws DimensionMismatch Latitude_LRF_State(
        B_by_lat[:, :, 1:2], chi_reference, q0, alpha, nλ)
    @test_throws ArgumentError Latitude_LRF_State(
        fill(NaN, nd, nd, nθ), chi_reference, q0, alpha, nλ)

    mktempdir() do directory
        filepath = joinpath(directory, "latitude_lrf_test.jld2")
        latitude = [-30.0, 0.0, 30.0]
        JLD2.jldsave(
            filepath;
            scheme = "regularized_log_latitude_v1",
            B_by_lat,
            chi_reference,
            q0_kg_kg = q0,
            alpha,
            latitude,
        )
        loaded = Load_LRF_State(filepath, nλ, nθ, nd; latitude = latitude)
        @test loaded isa Latitude_LRF_State
        @test loaded.B_by_lat == B_by_lat
        @test loaded.chi_reference == chi_reference
        @test loaded.q0 == q0
        @test loaded.alpha == alpha
        @test_throws ArgumentError Load_LRF_State(
            filepath, nλ, nθ, nd; latitude = reverse(latitude))
        @test_throws DimensionMismatch Load_LRF_State(filepath, nλ, nθ + 1, nd)
    end
end

@testset "LRF persistent storage and output registration" begin
    dyn_data = Dyn_Data("lrf_test", 1, 2, 4, 2, 2)
    @test size(dyn_data.grid_lrf_tendency) == (4, 2, 2)
    @test all(iszero, dyn_data.grid_lrf_tendency)

    live_data =
        JGCM.Variable_Mappings_Module.Get_Dyn_Var_Map(dyn_data, Val(:PrimitiveEquation))
    @test live_data[:lrf_dt] === dyn_data.grid_lrf_tendency

    metadata = JGCM.Output_Mappings_Module.Get_Var_Info(Val(:PrimitiveEquation))[:lrf_dt]
    @test metadata.nc_name == "lrf_dta_dt"
    @test metadata.units == "K s-1"
    @test metadata.dims == 3
end
