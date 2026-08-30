const SATURATION_FREEZE_TEMPERATURE = 273.16
const BM_MIN_PARCEL_TEMPERATURE = 173.16

"""
    Saturation_Vapor_Pressure(temperature)

Smithsonian saturation vapor pressure in Pa. Values are over ice below
253.16 K, over liquid above 273.16 K, and linearly blended in temperature
between those limits.
"""
function Saturation_Vapor_Pressure(temperature::Real)
    pressure, _ = _saturation_vapor_pressure_and_derivative(temperature)
    return pressure
end

function _saturation_vapor_pressure_and_derivative(temperature::Real)
    temperature = Float64(temperature)
    isfinite(temperature) && temperature > 0 ||
        throw(ArgumentError("temperature must be positive and finite"))

    t_freeze = SATURATION_FREEZE_TEMPERATURE
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

@inline function _validate_saturation_pressure(
    pressure::Float64,
    vapor_pressure::Float64,
    epsilon::Float64,
)
    isfinite(pressure) && pressure > 0.0 ||
        throw(ArgumentError("pressure must be positive and finite"))
    isfinite(epsilon) && 0.0 < epsilon < 1.0 ||
        throw(ArgumentError("epsilon must be finite and lie in (0, 1)"))
    vapor_pressure < pressure || throw(
        DomainError(
            (pressure = pressure, saturation_vapor_pressure = vapor_pressure),
            "saturation vapor pressure must be smaller than total pressure",
        ),
    )
    return nothing
end

"""
    Saturation_Mixing_Ratio(temperature, pressure, epsilon)

Exact saturation water-vapor mixing ratio (vapor mass per unit dry-air mass),
`epsilon * e_s / (pressure - e_s)`.
"""
@inline function Saturation_Mixing_Ratio(
    temperature::Real,
    pressure::Real,
    epsilon::Real,
)
    mixing_ratio, _ =
        _saturation_mixing_ratio_and_derivative(temperature, pressure, epsilon)
    return mixing_ratio
end

@inline function _saturation_mixing_ratio_and_derivative(
    temperature::Real,
    pressure::Real,
    epsilon::Real,
)
    pressure = Float64(pressure)
    epsilon = Float64(epsilon)
    es, des_dt = _saturation_vapor_pressure_and_derivative(temperature)
    _validate_saturation_pressure(pressure, es, epsilon)
    denominator = pressure - es
    mixing_ratio = epsilon * es / denominator
    derivative = epsilon * pressure * des_dt / denominator^2
    return mixing_ratio, derivative
end

"""
    Saturation_Specific_Humidity(temperature, pressure, epsilon)

Exact saturation specific humidity (vapor mass per unit moist-air mass),
`epsilon * e_s / (pressure - (1 - epsilon) * e_s)`.
"""
@inline function Saturation_Specific_Humidity(
    temperature::Real,
    pressure::Real,
    epsilon::Real,
)
    specific_humidity, _ =
        _saturation_specific_humidity_and_derivative(temperature, pressure, epsilon)
    return specific_humidity
end

@inline function _saturation_specific_humidity_and_derivative(
    temperature::Real,
    pressure::Real,
    epsilon::Real,
)
    pressure = Float64(pressure)
    epsilon = Float64(epsilon)
    es, des_dt = _saturation_vapor_pressure_and_derivative(temperature)
    _validate_saturation_pressure(pressure, es, epsilon)
    denominator = pressure - (1.0 - epsilon) * es
    specific_humidity = epsilon * es / denominator
    derivative = epsilon * pressure * des_dt / denominator^2
    return specific_humidity, derivative
end
