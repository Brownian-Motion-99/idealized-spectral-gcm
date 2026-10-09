# julia test/validation/f1_compare_20261006.jl BASELINE REPAIRED small|t42 REPORT
using Serialization, TOML

function compare_states(a,b)
    Dict(string(key)=>maximum(abs,a[key]-b[key]) for key in keys(a))
end
function main()
    length(ARGS)==4 || error("Usage: BASELINE REPAIRED small|t42 REPORT")
    baseline,repaired,mode,report = ARGS
    old = TOML.parsefile(joinpath(baseline,"results.toml"))
    new = TOML.parsefile(joinpath(repaired,"results.toml"))
    old["threads"]==new["threads"] || error("Thread counts differ")
    results = Dict{String,Any}("mode"=>mode,"threads"=>new["threads"])
    cases = mode=="small" ? ("hs","none") : ("bm_off","bm_on")
    for case in cases
        old1=deserialize(joinpath(baseline,"$(case)_step1.bin"))
        new1=deserialize(joinpath(repaired,"$(case)_step1.bin"))
        identical=all(old1[key]==new1[key] for key in keys(old1))
        # Independent FFTW.PATIENT plans can differ in operation order. The
        # physics-disabled case quantifies that cross-process roundoff floor.
        equivalent=all(isapprox(old1[key],new1[key];rtol=5e-13,atol=1e-18)
            for key in keys(old1))
        equivalent || error("First accepted prognostics changed beyond roundoff: $case")
        new[case]["invariant_failures"]==0 || error("Patched invariant failed: $case")
        results[case] = Dict("first_step_prognostics_identical"=>identical,
            "first_step_prognostics_equal_to_roundoff"=>equivalent,
            "first_step_max_abs_differences"=>compare_states(old1,new1),
            "baseline"=>old[case],"repaired"=>new[case])
        if mode=="small"
            old2=deserialize(joinpath(baseline,"$(case)_step2.bin"))
            new2=deserialize(joinpath(repaired,"$(case)_step2.bin"))
            results[case]["second_step_baseline_differences"] = compare_states(old2,new2)
            if case=="hs"
                ref2=deserialize(joinpath(baseline,"hs_reference_step2.bin"))
                matched=all(isapprox(new2[key],ref2[key];rtol=5e-13,atol=1e-18)
                    for key in keys(new2))
                matched || error("Patched trajectory differs from refreshed baseline")
                results[case]["second_step_matches_refreshed_reference"] = matched
                results[case]["second_step_reference_differences"] = compare_states(ref2,new2)
                maximum(abs,old2[:grid_u_c]-ref2[:grid_u_c])>1e-10 ||
                    error("Fixture did not exercise next-step consumption")
            else
                all(isapprox(old2[key],new2[key];rtol=5e-13,atol=1e-18)
                    for key in keys(old2)) || error("Physics-disabled trajectory changed")
                results[case]["second_step_prognostics_equal_to_roundoff"] = true
            end
        else
            new[case]["steps"]==144 || error("T42 run duration differs")
            case=="bm_on" && new[case]["bm_active_steps"]==0 && error("BM was inactive")
            new[case]["refresh_allocations_bytes"] < 128*64*20*8 ||
                error("Refresh allocated a full-grid-sized buffer")
            results[case]["median_step_time_ratio"] =
                new[case]["warmed_median_step_seconds"]/old[case]["warmed_median_step_seconds"]
        end
    end
    open(report,"w") do io; TOML.print(io,results); end
    println("$mode comparison passed; report: $report")
end
main()
