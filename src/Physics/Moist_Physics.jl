"""
    Moist_Physics!(...)

Sequentially apply Betts-Miller convection and large-scale condensation to a
private physics working state. Betts-Miller diagnostic rates are applied once
over one bounded physics substep; large-scale condensation then sees that
updated state directly.
"""
function Moist_Physics!(
    do_betts_miller::Bool,
    do_lscale_cond::Bool,
    bm_state::Union{Nothing,Betts_Miller_State},
    heating_fraction,
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

    total_precipitation .+= bm_precipitation

    if do_betts_miller
        @. temperature += effective_dt * bm_temperature_tendency
        @. humidity += effective_dt * bm_humidity_tendency
    end

    if do_lscale_cond
        Lscale_Cond!(
            atmo_data,
            temperature,
            humidity,
            p_full,
            p_half,
            effective_dt,
            heating_fraction,
            lscale_temperature_tendency,
            lscale_humidity_tendency,
            liquid_water_content,
            total_precipitation,
        )
    else
        fill!(lscale_temperature_tendency, 0.0)
        fill!(lscale_humidity_tendency, 0.0)
        fill!(liquid_water_content, 0.0)
    end
    return nothing
end
