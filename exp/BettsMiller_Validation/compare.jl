# Compare two completed runs of validate.jl; no additional package dependencies.
using Serialization, TOML, Printf

length(ARGS) == 3 || error("Usage: compare.jl RUN_600 RUN_300 OUTPUT_TOML")
coarse_dir, fine_dir, destination = abspath.(ARGS)
coarse = TOML.parsefile(joinpath(coarse_dir,"result.toml"))
fine = TOML.parsefile(joinpath(fine_dir,"result.toml"))
coarse["passed"] && fine["passed"] || error("Both runs must pass")
coarse["initial_state_sha256"] == fine["initial_state_sha256"] || error("Initial states differ")
coarse["days"] == fine["days"] && coarse["dt"] == 2fine["dt"] || error("Durations/timestep ratio differ")
coarse["source_sha256"] == fine["source_sha256"] || error("Runs used different source files")
c, f = deserialize(joinpath(coarse_dir,"final.bin")), deserialize(joinpath(fine_dir,"final.bin"))
mass_weight = 0.5 .* (diff(c.ph;dims=3) + diff(f.ph;dims=3)) .* reshape(c.weights,1,:,1)
weights = sum(mass_weight)
comparison = Dict{String,Any}()
for key in (:t,:q,:u,:v)
    difference = getfield(c,key) - getfield(f,key)
    comparison["$(key)_mass_weighted_rms"] = sqrt(sum(difference.^2 .* mass_weight) / weights)
    comparison["$(key)_max_absolute_difference"] = maximum(abs,difference)
end
ps_difference = c.ps-f.ps
comparison["ps_area_weighted_rms"] = sqrt(sum(ps_difference.^2 .* reshape(c.weights,1,:,1)) / (2size(c.ps,1)))
comparison["ps_max_absolute_difference"] = maximum(abs,ps_difference)
for name in ("bm_rain","lscale_rain","evaporation")
    key = "accumulated_$(name)_mm"
    comparison["$(name)_relative_difference"] = (coarse[key]-fine[key]) / max(abs(fine[key]),eps(Float64))
end
for report in (coarse,fine)
    s = report["statistics"]
    s["deep_samples"] + s["shallow_samples"] > 0 || error("Run never exercised BM adjustment")
    report["accumulated_lscale_rain_mm"] > 0 || error("Run never exercised condensation")
    report["accumulated_evaporation_mm"] > 0 || error("Run never exercised evaporation")
end
result = Dict("coarse"=>coarse,"fine"=>fine,"comparison"=>comparison,
    "identical_initial_states"=>true,"matched_source_files"=>true,
    "ledger_quantities"=>["water_kg_m2","cpT_plus_Lvq_plus_K_J_m2","dry_mass_kg_m2"],
    "ledger_method"=>"Previous-level leapfrog recurrence; process changes observed by replay matching production physics each step.",
    "energy_scope"=>"The ledger measures cp*T+Lv*q+K and records all observed changes. It is not a full moist-energy theorem; dynamics/fixers are grouped, and pressure/PBL terms are measured separately.",
    "passed"=>true)
mkpath(dirname(destination))
open(destination,"w") do io
    TOML.print(io,result;sorted=true)
end
println("Identical initial states and matched source files: passed")
for name in sort(collect(keys(comparison)))
    @printf("%s = %.8g\n",name,comparison[name])
end
