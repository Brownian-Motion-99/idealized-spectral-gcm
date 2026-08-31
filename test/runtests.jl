using Test
using JGCM

# 1. Math Kernel Tests
@testset "Math Kernels" begin
    include("test_Gauss_And_Legendre.jl")
    include("test_Spectral_Spherical_Mesh.jl")
    include("test_Grid_Tracer_Transport.jl")
    include("test_Transform_And_Helmholtz_Batching.jl")
    include("test_Diagnostic_Allocations.jl")
    include("test_Memory_And_Restart.jl")
    include("test_Initial_Conditions.jl")
    include("test_Lower_Boundary_Temperature.jl")
    include("test_LRF.jl")
    include("test_Betts_Miller.jl")
    include("test_Lscale_Cond.jl")
end

# 2. Driver Utility Tests
@testset "Driver Utilities" begin
    include("test_Driver.jl")
end

# 3. Physics parameterization tests
@testset "Physics Parameterizations" begin
    include("test_HS_Forcing.jl")
    include("test_PBL.jl")
    include("test_Virtual_Temperature_Dynamics.jl")
    include("test_Dry_Air_Adjustment.jl")
    include("test_Corrections.jl")
end

# 4. Integration Tests (Optional: Run short simulations)
@testset "Integration Tests" begin
    include("test_Coupling.jl")
end
