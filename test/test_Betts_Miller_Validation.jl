using Serialization
include(joinpath(@__DIR__, "support", "betts_miller_cases.jl"))

@testset "Betts-Miller seeded physical invariants" begin
    cases = BMValidation.ensemble()
    nd = length(first(cases).p_full)
    for mode in (:isca, :timescale), virtual in (false, true)
        state = Betts_Miller_State(nd; energy_correction = mode)
        atmo = BMValidation.atmosphere(nd; virtual)
        regimes = Set{Symbol}()
        for case in cases
            @testset "$(case.name) $mode virtual=$virtual" begin
                temperature, humidity = copy(case.temperature), copy(case.humidity)
                result = Betts_Miller_Column(state, atmo, temperature, humidity, case.p_full, case.p_half)
                push!(regimes, result.regime)
                @test temperature == case.temperature && humidity == case.humidity
                @test all(isfinite, result.temperature_tendency)
                @test all(isfinite, result.humidity_tendency)
                @test all(isfinite, result.parcel_temperature)
                @test isfinite(result.cape) && isfinite(result.cin)
                @test isfinite(result.precipitation) && result.precipitation >= 0
                @test result.active == (result.regime != :none)
                mass = diff(case.p_half) ./ atmo.grav
                water = result.humidity_tendency .* mass
                heat = atmo.cp_air .* result.temperature_tendency .* mass
                latent = atmo.Lv .* water
                tolerance = 64nd * eps(Float64)
                @test abs(result.precipitation + sum(water)) <=
                      tolerance * (result.precipitation + sum(abs.(water)))
                @test abs(sum(heat + latent)) <= tolerance * sum(abs.(heat) + abs.(latent))
                if result.active
                    @test 1 <= result.lzb <= result.adjustment_top <= nd
                    @test 0 < result.top_fraction <= 1
                    @test all(iszero, result.temperature_tendency[1:result.adjustment_top-1])
                    @test all(iszero, result.humidity_tendency[1:result.adjustment_top-1])
                    if result.regime == :shallow
                        @test result.precipitation == 0
                        @test abs(sum(water)) <= tolerance * sum(abs.(water))
                        @test abs(sum(heat)) <= tolerance * sum(abs.(heat) + abs.(latent))
                        @test result.reference_temperature ≈ temperature .+ state.tau .* result.temperature_tendency
                        @test result.reference_humidity ≈ humidity .+ state.tau .* result.humidity_tendency
                    else
                        @test result.adjustment_top == result.lzb && result.top_fraction == 1
                    end
                else
                    @test result.adjustment_top == 0 && result.top_fraction == 0
                    @test all(iszero, result.temperature_tendency)
                    @test all(iszero, result.humidity_tendency)
                end
                for dt in (600.0, state.tau / 2, state.tau)
                    final_q = humidity .+ dt .* result.humidity_tendency
                    final_t = temperature .+ dt .* result.temperature_tendency
                    @test all(isfinite, final_q) && all(0 .<= final_q .< 1)
                    @test all(isfinite, final_t) && all(final_t .> 0)
                    # Finite-step budgets use the unchanged column masses.
                    water_change = (final_q - humidity) .* mass
                    @test abs(sum(water_change) + dt * result.precipitation) <=
                          tolerance * (sum(abs.(humidity .* mass)) + dt * result.precipitation)
                end
                if startswith(case.name, "dry_neutral")
                    @test result.regime == :none && result.lcl == 0
                    @test result.parcel_temperature ≈ temperature rtol = 1.0e-14
                    @test result.parcel_mixing_ratio ≈ humidity ./ (1 .- humidity) rtol = 1.0e-14
                elseif case.name == "buoyant_to_finite_top"
                    @test result.lzb == 1 && result.cape > 0 && result.active
                elseif case.kind == :stable
                    # Thermal stability prevents adjustment; a dry environment
                    # can still give small positive virtual-temperature CAPE.
                    @test result.regime == :none
                    if result.cape > 0
                        top = result.lzb
                        preliminary_heat = (result.parcel_temperature[top:end] - temperature[top:end]) .* mass[top:end]
                        @test sum(preliminary_heat) <= tolerance * sum(abs.(preliminary_heat))
                    end
                end
            end
        end
        @test regimes == Set((:none, :deep, :shallow))
    end
end

function bm_kernel_allocations(state, atmo, case)
    work = state.work[1]
    args = (work, state, case.temperature, case.humidity, case.p_full, case.p_half,
        atmo.rdgas, atmo.rvgas, atmo.cp_air, atmo.Lv, atmo.grav, atmo.kappa, atmo.use_virtual_temperature)
    JGCM.Atmos_Param_Module._betts_miller_column!(args...)
    return @allocated JGCM.Atmos_Param_Module._betts_miller_column!(args...)
end

@testset "Betts-Miller LCL inference and allocation regression" begin
    nd = 30
    atmo = BMValidation.atmosphere(nd)
    state = Betts_Miller_State(nd)
    for case in (BMValidation.fixture("deep"), BMValidation.fixture("shallow"))
        theta = case.temperature[end] * (1.0e5 / case.p_full[end])^atmo.kappa
        r = case.humidity[end] / (1 - case.humidity[end])
        result = @inferred JGCM.Atmos_Param_Module._bm_lcl(
            theta, r, case.p_full[1], case.p_full[end], atmo.rdgas / atmo.rvgas, atmo.kappa,
        )
        @test result.found
        @test bm_kernel_allocations(state, atmo, case) <= 1024
    end
end

@testset "Betts-Miller fixed independent Isca ensemble" begin
    path = joinpath(@__DIR__, "fixtures", "betts_miller_isca_ensemble.tsv")
    rows = filter(line -> !startswith(line, "#"), readlines(path))
    @test length(rows) == 87 * 30
    for group in Iterators.partition(rows, 30)
        name = first(split(first(group)))
        data = permutedims(reduce(hcat, [parse.(Float64, split(line)[3:end]) for line in group]))
        pf, ph = data[:, 1], vcat(data[:, 2], data[end, 3])
        temperature, humidity = data[:, 4], data[:, 5]
        atmo = BMValidation.atmosphere(30)
        result = Betts_Miller_Column(Betts_Miller_State(30), atmo, temperature, humidity, pf, ph)
        @testset "$name" begin
            factor = 128 * 30 * eps(Float64)
            t_tol = factor * maximum(temperature) / BMValidation.TAU
            q_tol = factor * maximum(humidity) / BMValidation.TAU
            rain_tol = q_tol * sum(diff(ph) ./ atmo.grav)
            cape_tol = factor * atmo.rdgas * maximum(temperature) * sum(log.(ph[2:end] ./ ph[1:end-1]))
            # Strict exact-LCL control, which changes only the Isca LCL lookup.
            @test all(abs.(result.temperature_tendency - data[:, 13]) .<= t_tol)
            @test all(abs.(result.humidity_tendency - data[:, 14]) .<= q_tol)
            @test abs(result.precipitation - data[1, 15]) <= rain_tol
            @test abs(result.cape - data[1, 16]) <= cape_tol
            @test abs(result.cin - data[1, 17]) <= cape_tol
            @test result.lzb == Int(data[1, 19])
            control_lcl = Int(data[1, 18])
            if result.lcl != control_lcl
                r = humidity[end] / (1 - humidity[end])
                rs = Saturation_Mixing_Ratio(temperature[end], pf[end], atmo.rdgas / atmo.rvgas)
                @test abs(r - rs) <= 128eps(Float64) * abs(r)
                @test max(result.lcl, control_lcl) == 30 && abs(result.lcl - control_lcl) == 1
            end
            # Error against the original lookup is bounded by the measured
            # difference between the two independent Fortran references.
            @test all(abs.(result.temperature_tendency - data[:, 6]) .<= abs.(data[:, 13] - data[:, 6]) .+ t_tol)
            @test all(abs.(result.humidity_tendency - data[:, 7]) .<= abs.(data[:, 14] - data[:, 7]) .+ q_tol)
            @test abs(result.precipitation - data[1, 8]) <= abs(data[1, 15] - data[1, 8]) + rain_tol
            @test abs(result.cape - data[1, 9]) <= abs(data[1, 16] - data[1, 9]) + cape_tol
            @test abs(result.cin - data[1, 10]) <= abs(data[1, 17] - data[1, 10]) + cape_tol
        end
    end
end

@testset "Betts-Miller mixed-column grid and work reuse" begin
    cases = BMValidation.ensemble()
    temperature, humidity, pf, ph = BMValidation.grid_fields(cases)
    nlon, nlat, nd = size(temperature)
    for mode in (:isca, :timescale), virtual in (false, true)
        state = Betts_Miller_State(nd; energy_correction = mode)
        atmo = BMValidation.atmosphere(nd; virtual)
        dt, dq, rain = zeros(nlon, nlat, nd), zeros(nlon, nlat, nd), zeros(nlon, nlat, 1)
        buffers = [tuple((getfield(work, key) for key in fieldnames(typeof(work)))...) for work in state.work]
        for repetition in 1:2
            Betts_Miller!(state, atmo, temperature, humidity, pf, ph, dt, dq, rain)
            for j = 1:nlat, i = 1:nlon
                case = cases[mod1(i + (j - 1) * nlon, length(cases))]
                result = Betts_Miller_Column(state, atmo, case.temperature, case.humidity, case.p_full, case.p_half)
                @test dt[i, j, :] == result.temperature_tendency
                @test dq[i, j, :] == result.humidity_tendency
                @test rain[i, j, 1] == result.precipitation
            end
        end
        @test all(all(getfield(work, key) === buffer for (key, buffer) in zip(fieldnames(typeof(work)), original))
            for (work, original) in zip(state.work, buffers))
    end
end

@testset "Betts-Miller separate-process thread agreement" begin
    script = joinpath(@__DIR__, "validation", "thread_probe.jl")
    project = dirname(@__DIR__)
    mktempdir() do folder
        outputs = []
        for count in (1, 4)
            path = joinpath(folder, "threads_$count.bin")
            run(`$(Base.julia_cmd()) --startup-file=no --project=$project --threads=$count $script $path`)
            push!(outputs, deserialize(path))
        end
        @test outputs[1] == outputs[2]
    end
end
