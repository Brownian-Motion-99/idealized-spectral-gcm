using Test
using Random
using JGCM

const validate_physics_state! =
    JGCM.Atmos_Param_Module._validate_and_clean_physics_state!

@testset "Fused physics-state validation" begin
    rng = MersenneTwister(2028)
    u = randn(rng, 4, 3, 2)
    v = randn(rng, 4, 3, 2)
    temperature = 200.0 .+ 100.0 .* rand(rng, 4, 3, 2)
    humidity = 0.02 .* rand(rng, 4, 3, 2)
    humidity[1] = -1.0e-15
    expected_humidity = max.(humidity, 0.0)

    @test isnothing(validate_physics_state!(u, v, temperature, humidity))
    @test humidity == expected_humidity
    @test @allocated(validate_physics_state!(u, v, temperature, humidity)) == 0

    nonfinite_u = copy(u)
    nonfinite_u[1] = NaN
    @test_throws ErrorException validate_physics_state!(
        nonfinite_u, v, temperature, copy(expected_humidity),
    )

    invalid_temperature = copy(temperature)
    invalid_temperature[1] = 0.0
    @test_throws ErrorException validate_physics_state!(
        u, v, invalid_temperature, copy(expected_humidity),
    )

    negative_humidity = copy(expected_humidity)
    negative_humidity[1] = -1.0e-10
    @test_throws ErrorException validate_physics_state!(
        u, v, temperature, negative_humidity,
    )

    excessive_humidity = copy(expected_humidity)
    excessive_humidity[1] = 1.0
    @test_throws ErrorException validate_physics_state!(
        u, v, temperature, excessive_humidity,
    )
end
