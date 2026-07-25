"""
    Moist_Physics!(...)

Coordinate Betts-Miller convection and large-scale condensation. 
Betts-Miller is diagnosed first. Large-scale condensation then sees the 
temporary state obtained by applying the Betts-Miller rates over the 
effective Euler/leapfrog interval. Both schemes finally contribute additive
rates to the shared model tendencies.
"""
function Moist_Physics!(
    do_betts_miller::Bool,
    do_lscale_cond::Bool,
    bm_state::Union{Nothing,Betts_Miller_State},
    L,
    atmo_data::Atmo_Data,
    temperature::Array{Float64,3},
    humidity::Array{Float64,3},
    p_full::Array{Float64,3},
    p_half::Array{Float64,3},
    effective_dt::Real,
    bm_temperature_tendency::Array{Float64,3},
    bm_humidity_tendency::Array{Float64,3},
    bm_precipitation::Array{Float64,3},
    lscale_temperature_tendency::Array{Float64,3},
    lscale_humidity_tendency::Array{Float64,3},
    liquid_water_content::Array{Float64,3},
    total_temperature_tendency::Array{Float64,3},
    total_humidity_tendency::Array{Float64,3},
    total_precipitation::Array{Float64,3},
)
    if do_betts_miller
        isnothing(bm_state) && error("BM_state was not initialized before time integration")
        Betts_Miller!(
            bm_state,
            atmo_data,
            temperature,
            humidity,
            p_full,
            p_half,
            bm_temperature_tendency,
            bm_humidity_tendency,
            bm_precipitation,
        )
    else
        fill!(bm_temperature_tendency, 0.0)
        fill!(bm_humidity_tendency, 0.0)
        fill!(bm_precipitation, 0.0)
    end

    total_temperature_tendency .+= bm_temperature_tendency
    total_humidity_tendency .+= bm_humidity_tendency
    total_precipitation .+= bm_precipitation

    if do_lscale_cond
        Lscale_Cond!(
            atmo_data,
            temperature,
            humidity,
            p_full,
            p_half,
            effective_dt,
            L,
            bm_temperature_tendency,
            bm_humidity_tendency,
            lscale_temperature_tendency,
            lscale_humidity_tendency,
            liquid_water_content,
            total_precipitation,
        )
        total_temperature_tendency .+= lscale_temperature_tendency
        total_humidity_tendency .+= lscale_humidity_tendency
    else
        fill!(lscale_temperature_tendency, 0.0)
        fill!(lscale_humidity_tendency, 0.0)
        fill!(liquid_water_content, 0.0)
    end
    return nothing
end
