@testset "Driver timestep accounting" begin
    time_steps = JGCM.Driver.time_steps

    @test time_steps(86_400, 600) == 144
    @test time_steps(43_200, 600) == 72

    @test_throws ArgumentError time_steps(86_400, 0)
    @test_throws ArgumentError time_steps(0, 600)
    @test_throws ArgumentError time_steps(-600, 600)
    @test_throws ArgumentError time_steps(1_000, 600)
end

@testset "Driver progress metrics" begin
    halfway = JGCM.Driver.progress_metrics(72, 144, 72.0)
    @test halfway.completed_steps == 72
    @test halfway.total_steps == 144
    @test halfway.segment_progress == 0.5
    @test halfway.eta_seconds == 72.0

    not_started = JGCM.Driver.progress_metrics(0, 144, 0.0)
    @test not_started.segment_progress == 0.0
    @test isnothing(not_started.eta_seconds)

    complete = JGCM.Driver.progress_metrics(144, 144, 120.0)
    @test complete.segment_progress == 1.0
    @test complete.eta_seconds == 0.0

    @test_throws ArgumentError JGCM.Driver.progress_metrics(0, 0, 0.0)
    @test_throws ArgumentError JGCM.Driver.progress_metrics(145, 144, 120.0)
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
