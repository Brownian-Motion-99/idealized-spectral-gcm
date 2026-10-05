using JGCM, Printf
include(joinpath(@__DIR__, "..", "support", "betts_miller_cases.jl"))

output = only(ARGS)
mkpath(output)
cases = BMValidation.ensemble()
open(joinpath(output, "columns.txt"), "w") do io
    println(io, length(cases), " ", length(first(cases).p_full))
    for case in cases
        println(io, case.reference ? 1 : 0)
        for k in eachindex(case.p_full)
            @printf(io, "%.17g %.17g %.17g %.17g %.17g\n", case.p_full[k], case.p_half[k],
                case.p_half[k+1], case.temperature[k], case.humidity[k])
        end
    end
end
open(joinpath(output, "julia.txt"), "w") do io
    state = Betts_Miller_State(length(first(cases).p_full))
    atmo = BMValidation.atmosphere(state.nd)
    for case in cases
        result = Betts_Miller_Column(state, atmo, case.temperature, case.humidity, case.p_full, case.p_half)
        println(io, case.name, " ", result.regime, " ", result.lcl, " ", result.lzb,
            " ", result.adjustment_top, " ", result.top_fraction, " ", result.cape,
            " ", result.cin, " ", result.precipitation)
        for k in eachindex(case.p_full)
            @printf(io, "%.17g %.17g %.17g %.17g %.17g\n", result.temperature_tendency[k],
                result.humidity_tendency[k], result.reference_temperature[k],
                result.reference_humidity[k], result.parcel_temperature[k])
        end
    end
end
