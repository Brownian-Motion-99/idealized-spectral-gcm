# Run unchanged with either the pre-fix or corrected JGCM project to compare
# fixed inputs. Outputs are measurements, not test expectations.
using JGCM, Printf
include(joinpath(@__DIR__,"..","..","test","support","betts_miller_cases.jl"))
println("name\tactive\tcape_J_kg\train_mm_day\tsurface_Tdot_K_s\tmax_abs_qdot_s\tlcl\tlzb\tmax_parcel_error_K")
function emit(name, temperature, humidity, pf, ph; dry = false)
    atmo = BMValidation.atmosphere(length(pf))
    result = Betts_Miller_Column(Betts_Miller_State(length(pf)),atmo,temperature,humidity,pf,ph)
    @printf("%s\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%d\t%d\t%.17g\n",
        name,string(result.active),result.cape,86400result.precipitation,
        result.temperature_tendency[end],maximum(abs,result.humidity_tendency),
        result.lcl,result.lzb,dry ? maximum(abs.(result.parcel_temperature-temperature)) : 0.0)
end
for name in ("deep","shallow","virtual")
    c = BMValidation.fixture(name)
    emit(name,c.temperature,c.humidity,c.p_full,c.p_half)
end
ph = [10000.0,25000.0,40000.0,55000.0,70000.0,85000.0,100000.0]
pf = 0.5 .* (ph[1:end-1]+ph[2:end])
t = 300.0 .* (pf./pf[end]).^(2/7)
emit("dry_neutral",t,fill(1.0e-12,6),pf,ph;dry=true)
