using JGCM, Serialization
include(joinpath(@__DIR__, "..", "support", "betts_miller_cases.jl"))
cases = BMValidation.ensemble()
results = Dict((mode, virtual) => BMValidation.grid_outputs(cases; mode, virtual)
    for mode in (:isca, :timescale), virtual in (false, true))
serialize(only(ARGS), results)
