@testset "Driver timestep accounting" begin
    remaining_time_steps = JGCM.Driver.remaining_time_steps

    @test remaining_time_steps(0, 86_400, 600) == 144
    @test remaining_time_steps(43_200, 86_400, 600) == 72
    @test remaining_time_steps(86_400, 86_400, 600) == 0

    @test_throws ArgumentError remaining_time_steps(0, 86_400, 0)
    @test_throws ArgumentError remaining_time_steps(86_400, 43_200, 600)
    @test_throws ArgumentError remaining_time_steps(1, 1_000, 600)
end

@testset "Driver progress metrics" begin
    cold = JGCM.Driver.progress_metrics(43_200, 0, 86_400, 600, 72.0)
    @test cold.completed_steps == 72
    @test cold.total_steps == 144
    @test cold.segment_progress == 0.5
    @test cold.overall_progress == 0.5
    @test cold.eta_seconds == 72.0

    restart = JGCM.Driver.progress_metrics(64_800, 43_200, 86_400, 600, 30.0)
    @test restart.completed_steps == 36
    @test restart.total_steps == 72
    @test restart.segment_progress == 0.5
    @test restart.overall_progress == 0.75
    @test restart.eta_seconds == 30.0
end

@testset "Driver end-of-run metrics" begin
    metrics = JGCM.Driver.run_metrics(43_200, 86_400, 72, 36.0, 86_400)
    @test metrics.simulated_days == 0.5
    @test metrics.seconds_per_step == 0.5
    @test metrics.simulated_days_per_wall_day == 1200.0

    no_steps = JGCM.Driver.run_metrics(86_400, 86_400, 0, 0.0, 86_400)
    @test no_steps.simulated_days == 0.0
    @test isnothing(no_steps.seconds_per_step)
    @test isnothing(no_steps.simulated_days_per_wall_day)
end
