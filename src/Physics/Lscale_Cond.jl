using Base.Threads

@inline _lscale_heating_scale(scale::Real, i::Int, j::Int) = Float64(scale)
@inline _lscale_heating_scale(scale::AbstractArray{<:Real,2}, i::Int, j::Int) =
    Float64(scale[i, j])

function _validate_lscale_heating_scale(scale, nλ::Int, nθ::Int)
    if scale isa Real
        value = Float64(scale)
        isfinite(value) && 0.0 <= value <= 1.0 ||
            throw(ArgumentError("L must be finite and lie in [0, 1]"))
    elseif scale isa AbstractArray{<:Real,2}
        size(scale) == (nλ, nθ) ||
            throw(DimensionMismatch("array-valued L must have size (nλ, nθ)"))
        all(value -> isfinite(value) && 0.0 <= value <= 1.0, scale) ||
            throw(ArgumentError("every value of L must be finite and lie in [0, 1]"))
    else
        throw(ArgumentError("L must be a real scalar or a two-dimensional real array"))
    end
    return nothing
end

"""
    Lscale_Cond!(
        atmo_data, temperature, humidity, p_full, p_half, effective_dt, L,
        prior_temperature_tendency, prior_humidity_tendency,
        temperature_tendency, humidity_tendency,
        liquid_water_content, precipitation,
    )

Diagnose large-scale condensation after applying the supplied prior tendencies
over `effective_dt`. The prior tendencies are Betts-Miller tendencies. Large-scale
condensation outputs are rates, use the model-wide signed tendency convention
(`humidity_tendency < 0` for condensation), and are overwritten. Precipitation
is a positive flux and is added to the supplied column accumulator.

`L` scales latent heating only. It does not alter condensed water or
precipitation and intentionally does not appear in the saturation-adjustment
denominator.
"""
function Lscale_Cond!(
    atmo_data::Atmo_Data,
    temperature::Array{Float64,3},
    humidity::Array{Float64,3},
    p_full::Array{Float64,3},
    p_half::Array{Float64,3},
    effective_dt::Real,
    L,
    prior_temperature_tendency::Array{Float64,3},
    prior_humidity_tendency::Array{Float64,3},
    temperature_tendency::Array{Float64,3},
    humidity_tendency::Array{Float64,3},
    liquid_water_content::Array{Float64,3},
    precipitation::Array{Float64,3},
)
    size(temperature) == size(humidity) == size(p_full) ||
        throw(DimensionMismatch("large-scale-condensation full-level fields must match"))
    size(temperature) ==
    size(prior_temperature_tendency) ==
    size(prior_humidity_tendency) ==
    size(temperature_tendency) ==
    size(humidity_tendency) ==
    size(liquid_water_content) ||
        throw(DimensionMismatch("large-scale-condensation tendency fields must match"))

    nλ, nθ, nd = size(temperature)
    size(p_half) == (nλ, nθ, nd + 1) ||
        throw(DimensionMismatch("large-scale-condensation p_half has incorrect size"))
    size(precipitation) == (nλ, nθ, 1) || throw(
        DimensionMismatch("large-scale-condensation precipitation has incorrect size"),
    )
    (nλ, nθ, nd) == (atmo_data.nλ, atmo_data.nθ, atmo_data.nd) ||
        throw(DimensionMismatch("large-scale-condensation fields do not match atmosphere"))

    effective_dt = Float64(effective_dt)
    isfinite(effective_dt) && effective_dt > 0 ||
        throw(ArgumentError("effective_dt must be positive and finite"))
    _validate_lscale_heating_scale(L, nλ, nθ)

    cp = atmo_data.cp_air
    lv = atmo_data.Lv
    epsilon = atmo_data.rdgas / atmo_data.rvgas
    lv_over_cp = lv / cp
    grav = atmo_data.grav

    fill!(temperature_tendency, 0.0)
    fill!(humidity_tendency, 0.0)
    fill!(liquid_water_content, 0.0)

    # Each thread owns complete columns, so precipitation accumulation is
    # deterministic and free of the level-wise data race in the old kernel.
    @threads for j = 1:nθ
        for i = 1:nλ
            column_precipitation_amount = 0.0
            heating_scale = _lscale_heating_scale(L, i, j)

            for k = 1:nd
                t_star =
                    temperature[i, j, k] +
                    effective_dt * prior_temperature_tendency[i, j, k]
                q_star = humidity[i, j, k] + effective_dt * prior_humidity_tendency[i, j, k]
                pressure = p_full[i, j, k]

                isfinite(t_star) && t_star > 0 || throw(
                    ArgumentError(
                        "post-convection temperature must be positive and finite",
                    ),
                )
                isfinite(q_star) && 0.0 <= q_star < 1.0 || throw(
                    ArgumentError(
                        "post-convection specific humidity must satisfy 0 <= q < 1",
                    ),
                )

                q_sat, dq_sat_dt = _bm_saturation_specific_humidity_and_derivative(
                    t_star,
                    pressure,
                    epsilon,
                )

                if q_star > q_sat && q_sat > 0.0
                    humidity_increment = (q_sat - q_star) / (1.0 + lv_over_cp * dq_sat_dt)
                    temperature_increment = -heating_scale * lv_over_cp * humidity_increment

                    humidity_rate = humidity_increment / effective_dt
                    temperature_rate = temperature_increment / effective_dt
                    humidity_tendency[i, j, k] = humidity_rate
                    temperature_tendency[i, j, k] = temperature_rate
                    liquid_water_content[i, j, k] = -humidity_rate

                    layer_mass = (p_half[i, j, k+1] - p_half[i, j, k]) / grav
                    layer_mass >= 0 ||
                        throw(ArgumentError("p_half must increase from top to surface"))
                    column_precipitation_amount -= humidity_increment * layer_mass
                end
            end

            precipitation[i, j, 1] += column_precipitation_amount / effective_dt
        end
    end
    return nothing
end
