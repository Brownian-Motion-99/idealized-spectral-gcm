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
