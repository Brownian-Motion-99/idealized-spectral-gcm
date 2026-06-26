using Test
using JGCM

# 1. Math Kernel Tests
@testset "Math Kernels" begin
    include("test_Gauss_And_Legendre.jl")
    include("test_Spectral_Spherical_Mesh.jl")
end

# 2. Integration Tests (Optional: Run short simulations)
@testset "Integration Tests" begin
    # You could potentially include a very short Barotropic run here
    # to ensure the full pipeline works without crashing.
    # include("../exp/Barotropic/Barotropic.jl") 
end
