const BM_FREEZE_TEMPERATURE = 273.16
const BM_MIN_PARCEL_TEMPERATURE = 173.16

"""
    Betts_Miller_Saturation_Vapor_Pressure(temperature)

Smithsonian saturation vapor pressure in Pa. Values are over ice below
253.16 K, over liquid above 273.16 K, and linearly blended in temperature
between those limits.
"""
function Betts_Miller_Saturation_Vapor_Pressure(temperature::Real)
    pressure, _ = _bm_saturation_vapor_pressure_and_derivative(temperature)
    return pressure
end

function _bm_saturation_vapor_pressure_and_derivative(temperature::Real)
    temperature = Float64(temperature)
    isfinite(temperature) && temperature > 0 ||
        throw(ArgumentError("temperature must be positive and finite"))

    t_freeze = BM_FREEZE_TEMPERATURE
    t_water_base = t_freeze + 100.0
    log_ten = log(10.0)

    es_ice = 0.0
    des_ice = 0.0
    if temperature < t_freeze
        ratio = t_freeze / temperature
        x =
            -9.09718 * (ratio - 1.0) - 3.56654 * log10(ratio) +
            0.876793 * (1.0 - temperature / t_freeze) +
            log10(610.71)
        es_ice = 10.0^x
        dx_dt =
            9.09718 * t_freeze / temperature^2 + 3.56654 / (temperature * log_ten) -
            0.876793 / t_freeze
        des_ice = log_ten * es_ice * dx_dt
    end

    es_liquid = 0.0
    des_liquid = 0.0
    if temperature > t_freeze - 20.0
        ratio = t_water_base / temperature
        liquid_power_1 = (1.0 - temperature / t_water_base) * 11.344
        liquid_power_2 = (ratio - 1.0) * (-3.49149)
        power_1 = 10.0^liquid_power_1
        power_2 = 10.0^liquid_power_2
        x =
            -7.90298 * (ratio - 1.0) + 5.02808 * log10(ratio) -
            1.3816e-7 * (power_1 - 1.0) +
            8.1328e-3 * (power_2 - 1.0) +
            log10(101324.60)
        es_liquid = 10.0^x
        dx_dt =
            7.90298 * t_water_base / temperature^2 - 5.02808 / (temperature * log_ten) +
            1.3816e-7 * log_ten * 11.344 / t_water_base * power_1 +
            8.1328e-3 * log_ten * 3.49149 * t_water_base / temperature^2 * power_2
        des_liquid = log_ten * es_liquid * dx_dt
    end

    if temperature <= t_freeze - 20.0
        return es_ice, des_ice
    elseif temperature >= t_freeze
        return es_liquid, des_liquid
    end

    ice_weight = (t_freeze - temperature) / 20.0
    es = ice_weight * es_ice + (1.0 - ice_weight) * es_liquid
    des =
        ice_weight * des_ice + (1.0 - ice_weight) * des_liquid + (es_liquid - es_ice) / 20.0
    return es, des
end

@inline function _bm_saturation_mixing_ratio(
    temperature::Float64,
    pressure::Float64,
    epsilon::Float64,
)
    return epsilon * Betts_Miller_Saturation_Vapor_Pressure(temperature) / pressure
end

@inline function _bm_saturation_specific_humidity_and_derivative(
    temperature::Float64,
    pressure::Float64,
    epsilon::Float64,
)
    isfinite(pressure) && pressure > 0 ||
        throw(ArgumentError("pressure must be positive and finite"))
    es, des_dt = _bm_saturation_vapor_pressure_and_derivative(temperature)
    saturation_mixing_ratio = epsilon * es / pressure
    denominator = 1.0 + saturation_mixing_ratio
    saturation_specific_humidity = saturation_mixing_ratio / denominator
    derivative = epsilon * des_dt / pressure / denominator^2
    return saturation_specific_humidity, derivative
end
