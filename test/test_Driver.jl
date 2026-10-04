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

@testset "Driver Betts-Miller configuration" begin
    make_state = JGCM.Driver.betts_miller_state
    params = Dict{String,Any}("do_Betts_Miller" => true)
    default = make_state(params, 6, 600, true)
    @test default.energy_correction == :isca
    @test default.tau == 7200.0
    @test default.relative_humidity == 0.8
    for mode in (:isca, :timescale, "isca", "timescale")
        configured = merge(params, Dict(
            "bm_energy_correction" => mode,
            "bm_tau" => 1800.0,
            "bm_relative_humidity" => 0.7,
        ))
        state = make_state(configured, 6, 1800, true)
        @test state.energy_correction == Symbol(mode)
        @test state.tau == 1800.0
        @test state.relative_humidity == 0.7
    end
    for invalid in (:unknown, "ISCA", 1, nothing)
        @test_throws ArgumentError make_state(
            merge(params, Dict("bm_energy_correction" => invalid)), 6, 600, true,
        )
    end
    @test_throws ArgumentError make_state(params, 6, 600, false)
    @test_throws ArgumentError make_state(params, 6, 7201, true)
    @test_throws ArgumentError make_state(merge(params, Dict("bm_tau" => 0.0)), 6, 600, true)
    @test_throws ArgumentError make_state(merge(params, Dict("bm_relative_humidity" => 1.1)), 6, 600, true)
    @test isnothing(make_state(Dict{String,Any}(), 6, 600, false))
end
