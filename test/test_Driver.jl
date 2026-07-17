@testset "Driver timestep accounting" begin
    remaining_time_steps = JGCM.Driver.remaining_time_steps

    @test remaining_time_steps(0, 86_400, 600) == 144
    @test remaining_time_steps(43_200, 86_400, 600) == 72
    @test remaining_time_steps(86_400, 86_400, 600) == 0

    @test_throws ArgumentError remaining_time_steps(0, 86_400, 0)
    @test_throws ArgumentError remaining_time_steps(86_400, 43_200, 600)
    @test_throws ArgumentError remaining_time_steps(1, 1_000, 600)
end
