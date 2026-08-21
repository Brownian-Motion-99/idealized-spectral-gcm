using Test
using JGCM

# 1. Math Kernel Tests
@testset "Math Kernels" begin
    include("test_Gauss_And_Legendre.jl")
    include("test_Spectral_Spherical_Mesh.jl")
end

# 2. Driver Utility Tests
@testset "Driver Utilities" begin
    include("test_Driver.jl")
end

# 3. Physics parameterization tests
@testset "Physics Parameterizations" begin
    include("test_PBL.jl")
    include("test_Corrections.jl")
end

# 4. Integration Tests (Optional: Run short simulations)
@testset "Integration Tests" begin
    # You could potentially include a very short Barotropic run here
    # to ensure the full pipeline works without crashing.
    # include("../exp/Barotropic/Barotropic.jl") 
end
