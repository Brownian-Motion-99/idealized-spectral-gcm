using Base.Threads

mutable struct Betts_Miller_Work
    parcel_temperature::Vector{Float64}
    parcel_mixing_ratio::Vector{Float64}
    parcel_saturation_mixing_ratio::Vector{Float64}
    reference_temperature::Vector{Float64}
    reference_humidity::Vector{Float64}
    temperature_tendency::Vector{Float64}
    humidity_tendency::Vector{Float64}
    moisture_suffix::Vector{Float64}
end

function Betts_Miller_Work(nd::Int)
    nd > 1 || throw(ArgumentError("Betts-Miller requires at least two vertical levels"))
    return Betts_Miller_Work(
        zeros(nd),
        zeros(nd),
        zeros(nd),
        zeros(nd),
        zeros(nd),
        zeros(nd),
        zeros(nd),
        zeros(nd + 1),
    )
end

"""
    Betts_Miller_State(nd; tau=7200.0, relative_humidity=0.8, energy_correction=:isca)

Time-independent configuration and per-thread work storage for the Betts-Miller
convective adjustment. The scheme returns rates; the model time step is not part
of the column calculation. Deep-convection enthalpy closure uses the Isca
reference-temperature correction by default; `:timescale` retains the previous
scaling of the larger precipitation-equivalent adjustment.
"""
struct Betts_Miller_State
    tau::Float64
    relative_humidity::Float64
    energy_correction::Symbol
    nd::Int
    work::Vector{Betts_Miller_Work}
end

function Betts_Miller_State(
    nd::Int;
    tau::Real = 7200.0,
    relative_humidity::Real = 0.8,
    energy_correction = :isca,
)
    tau = Float64(tau)
    relative_humidity = Float64(relative_humidity)
    isfinite(tau) && tau > 0 || throw(ArgumentError("bm_tau must be positive and finite"))
    isfinite(relative_humidity) && 0 < relative_humidity <= 1 ||
        throw(ArgumentError("bm_relative_humidity must lie in (0, 1]"))
    (energy_correction isa Symbol || energy_correction isa AbstractString) ||
        throw(ArgumentError("bm_energy_correction must be :isca or :timescale"))
    energy_correction = Symbol(energy_correction)
    energy_correction in (:isca, :timescale) ||
        throw(ArgumentError("bm_energy_correction must be :isca or :timescale"))
    return Betts_Miller_State(
        tau,
        relative_humidity,
        energy_correction,
        nd,
        [Betts_Miller_Work(nd) for _ = 1:Threads.nthreads()],
    )
end

@inline _bm_layer_mass(p_half, k, grav) =
    (Float64(p_half[k+1]) - Float64(p_half[k])) / grav

function _bm_precipitation_integrals(tdot, qdot, p_half, top, surface, cp, lv, grav)
    moisture = 0.0
    thermal = 0.0
    moisture_scale = 0.0
    thermal_scale = 0.0
    for k = top:surface
        mass = _bm_layer_mass(p_half, k, grav)
        drying = -qdot[k] * mass
        heating = (cp / lv) * tdot[k] * mass
        moisture += drying
        thermal += heating
        moisture_scale += abs(drying)
        thermal_scale += abs(heating)
    end
    # Bound summation roundoff using the actual term magnitudes, without a
    # physical precipitation or CAPE threshold.
    relative_tolerance = 8 * (surface - top + 1) * eps(Float64)
    return (;
        moisture,
        thermal,
        moisture_tolerance = relative_tolerance * moisture_scale,
        thermal_tolerance = relative_tolerance * thermal_scale,
    )
end

function _bm_conserve_enthalpy!(work, p_half, top, surface, cp, lv, grav, tau)
    thermal = 0.0
    moisture = 0.0
    mass = 0.0
    for k = top:surface
        layer_mass = _bm_layer_mass(p_half, k, grav)
        thermal += work.temperature_tendency[k] * layer_mass
        moisture += work.humidity_tendency[k] * layer_mass
        mass += layer_mass
    end
    # Sum thermal and moisture terms separately to retain accuracy when large
    # opposing shallow moisture fluxes have a zero column integral.
    correction = -(thermal + (lv / cp) * moisture) / mass
    for k = top:surface
        work.temperature_tendency[k] += correction
        work.reference_temperature[k] += tau * correction
    end
    return correction
end

function _bm_deep_convection!(work, state, p_half, top, surface, cp, lv, grav, integrals)
    if integrals.moisture > integrals.thermal
        scale = integrals.thermal / integrals.moisture
        @views work.humidity_tendency[top:surface] .*= scale
        return integrals.thermal
    elseif state.energy_correction == :timescale
        scale = integrals.moisture / integrals.thermal
        @views work.temperature_tendency[top:surface] .*= scale
    else
        _bm_conserve_enthalpy!(work, p_half, top, surface, cp, lv, grav, state.tau)
    end
    return integrals.moisture
end

@inline function _bm_reference_humidity(saturation_mixing_ratio, pressure, rh, epsilon)
    vapor_pressure = rh * pressure * saturation_mixing_ratio / (epsilon + saturation_mixing_ratio)
    return epsilon * vapor_pressure / (pressure - (1.0 - epsilon) * vapor_pressure)
end

function _bm_reset_adjustment!(work, temperature, humidity, top, surface)
    for k = top:surface
        work.temperature_tendency[k] = 0.0
        work.humidity_tendency[k] = 0.0
        work.reference_temperature[k] = Float64(temperature[k])
        work.reference_humidity[k] = max(Float64(humidity[k]), 0.0)
    end
    return nothing
end

function _bm_shallow_convection!(
    work, state, temperature, humidity, p_half, lzb, surface, cp, lv, grav, integrals,
)
    top = lzb
    fraction = 1.0
    found = abs(integrals.moisture) <= integrals.moisture_tolerance
    if !found
        # Sum suffixes from the surface, so removing a large upper contribution
        # cannot erase the small lower-column integral used to solve for f.
        suffix = work.moisture_suffix
        suffix[surface+1] = 0.0
        for k = surface:-1:lzb
            suffix[k] =
                suffix[k+1] + work.humidity_tendency[k] * _bm_layer_mass(p_half, k, grav)
        end
        for k = lzb:surface
            lower = suffix[k+1]
            transition = work.humidity_tendency[k] * _bm_layer_mass(p_half, k, grav)
            if transition > 0 && lower <= 0 && suffix[k] >= 0
                fraction = clamp(-lower / transition, 0.0, 1.0)
                top = k
                if fraction == 0.0
                    # A zero fraction is an exact interface crossing. Exclude
                    # this cell from both the tendencies and correction mass.
                    top += 1
                    fraction = 1.0
                end
                found = true
                break
            end
        end
    end

    # A single cell has no nontrivial transport that conserves both water and
    # enthalpy. This also covers a root requiring a zero surface fraction.
    if !found || top >= surface
        _bm_reset_adjustment!(work, temperature, humidity, lzb, surface)
        return (; active = false, adjustment_top = 0, top_fraction = 0.0)
    end

    _bm_reset_adjustment!(work, temperature, humidity, lzb, top - 1)
    work.temperature_tendency[top] *= fraction
    work.humidity_tendency[top] *= fraction
    temperature_scale = 0.0
    for k = top:surface
        temperature_scale = max(temperature_scale, abs(work.temperature_tendency[k]))
        # These effective references include fractional penetration once.
        work.reference_temperature[k] =
            Float64(temperature[k]) + state.tau * work.temperature_tendency[k]
        work.reference_humidity[k] =
            max(Float64(humidity[k]), 0.0) + state.tau * work.humidity_tendency[k]
    end
    # Use full cell masses after f has been applied to the cell rates, as in
    # Isca's discrete closure. Do not apply f again to the correction mass.
    correction = _bm_conserve_enthalpy!(work, p_half, top, surface, cp, lv, grav, state.tau)
    temperature_tolerance =
        8 * (surface - top + 1) * eps(Float64) * (temperature_scale + abs(correction))
    active = any(k -> abs(work.temperature_tendency[k]) > temperature_tolerance, top:surface) ||
             any(!iszero, work.humidity_tendency)
    if !active
        # A uniform preliminary warming with no moisture transport can cancel
        # to roundoff. Do not report that numerical residue as convection.
        _bm_reset_adjustment!(work, temperature, humidity, top, surface)
        return (; active = false, adjustment_top = 0, top_fraction = 0.0)
    end
    return (; active, adjustment_top = top, top_fraction = fraction)
end

function _bm_close_adjustment!(work, state, temperature, humidity, p_half, top, surface, cp, lv, grav)
    integrals = _bm_precipitation_integrals(
        work.temperature_tendency, work.humidity_tendency, p_half, top, surface, cp, lv, grav,
    )
    if integrals.thermal <= integrals.thermal_tolerance
        _bm_reset_adjustment!(work, temperature, humidity, top, surface)
        return (;
            active = false, regime = :none, adjustment_top = 0, top_fraction = 0.0,
            precipitation = 0.0,
        )
    elseif integrals.moisture <= integrals.moisture_tolerance
        shallow = _bm_shallow_convection!(
            work, state, temperature, humidity, p_half, top, surface, cp, lv, grav, integrals,
        )
        return merge(shallow, (;
            regime = shallow.active ? :shallow : :none, precipitation = 0.0,
        ))
    end
    precipitation =
        _bm_deep_convection!(work, state, p_half, top, surface, cp, lv, grav, integrals)
    return (; active = true, regime = :deep, adjustment_top = top, top_fraction = 1.0, precipitation)
end

function _validate_bm_column(
    temperature::AbstractVector{<:Real},
    humidity::AbstractVector{<:Real},
    p_full::AbstractVector{<:Real},
    p_half::AbstractVector{<:Real},
)
    nd = length(p_full)
    length(temperature) == nd || throw(DimensionMismatch("temperature and p_full differ"))
    length(humidity) == nd || throw(DimensionMismatch("humidity and p_full differ"))
    length(p_half) == nd + 1 ||
        throw(DimensionMismatch("length(p_half) must equal length(p_full) + 1"))
    nd > 1 || throw(ArgumentError("Betts-Miller requires at least two vertical levels"))

    all(isfinite, temperature) || throw(ArgumentError("temperature must be finite"))
    all(q -> isfinite(q) && q < 1, humidity) || throw(
        ArgumentError(
            "specific humidity must be finite and less than 1; " *
            "column extrema are $(extrema(humidity))",
        ),
    )
    all(p -> isfinite(p) && p > 0, p_full) ||
        throw(ArgumentError("full-level pressure must be positive and finite"))
    all(p -> isfinite(p) && p >= 0, p_half) ||
        throw(ArgumentError("half-level pressure must be nonnegative and finite"))
    all(diff(p_full) .> 0) ||
        throw(ArgumentError("p_full must increase monotonically from top to surface"))
    all(diff(p_half) .> 0) ||
        throw(ArgumentError("p_half must increase monotonically from top to surface"))
    for k = 1:nd
        p_half[k] < p_full[k] < p_half[k+1] || throw(
            ArgumentError(
                "full-level pressure must lie between its surrounding interfaces",
            ),
        )
    end
    return nothing
end

@inline function _bm_layer_log_pressure(
    p_full::AbstractVector{<:Real},
    p_half::AbstractVector{<:Real},
    k::Int,
)
    # A zero-pressure model top makes the reference log(p₂/p₁) singular.
    # Use the uppermost full level as the finite top bound for that half layer.
    upper_pressure = p_half[k] > 0 ? p_half[k] : p_full[k]
    return log(p_half[k+1] / upper_pressure)
end

@inline function _bm_buoyancy(
    parcel_temperature::Float64,
    parcel_mixing_ratio::Float64,
    temperature::Real,
    humidity::Real,
    virtual_coefficient::Float64,
)
    parcel_humidity = parcel_mixing_ratio / (1.0 + parcel_mixing_ratio)
    parcel_virtual = parcel_temperature * (1.0 + virtual_coefficient * parcel_humidity)
    environment_virtual =
        Float64(temperature) * (1.0 + virtual_coefficient * max(Float64(humidity), 0.0))
    buoyancy = parcel_virtual - environment_virtual
    # A neutral parcel must not acquire an LFC from floating-point roundoff.
    tolerance = 8 * eps(Float64) * max(abs(parcel_virtual), abs(environment_virtual))
    return abs(buoyancy) <= tolerance ? 0.0 : buoyancy
end

function _bm_lcl(
    theta0::Float64,
    r0::Float64,
    p_top::Float64,
    p_surface::Float64,
    epsilon::Float64,
    kappa::Float64,
)
    pstar = 1.0e5
    residual(logp) = begin
        pressure = exp(logp)
        parcel_temperature = theta0 * (pressure / pstar)^kappa
        Saturation_Mixing_Ratio(parcel_temperature, pressure, epsilon) - r0
    end

    lo = log(p_top)
    hi = log(p_surface)
    temperature_top = theta0 * (p_top / pstar)^kappa
    saturation_top = Saturation_Mixing_Ratio(temperature_top, p_top, epsilon)
    tolerance = 32 * eps(Float64) * max(saturation_top, r0)
    if r0 <= 0 || saturation_top - r0 > tolerance
        return (; found = false, pressure = p_top, temperature = temperature_top)
    elseif abs(saturation_top - r0) <= tolerance
        return (; found = true, pressure = p_top, temperature = temperature_top)
    end

    for _ = 1:80
        mid = 0.5 * (lo + hi)
        if residual(mid) > 0
            hi = mid
        else
            lo = mid
        end
    end
    pressure = exp(0.5 * (lo + hi))
    return (; found = true, pressure, temperature = theta0 * (pressure / pstar)^kappa)
end

@inline function _bm_moist_derivative(
    temperature::Float64,
    mixing_ratio::Float64,
    kappa::Float64,
    cp::Float64,
    lv::Float64,
    rv::Float64,
)
    numerator = kappa * temperature + (lv / cp) * mixing_ratio
    denominator = 1.0 + lv^2 * mixing_ratio / (cp * rv * temperature^2)
    return numerator / denominator
end

function _bm_moist_rk2(
    temperature_a::Float64,
    mixing_ratio_a::Float64,
    pressure_a::Float64,
    pressure_b::Float64,
    epsilon::Float64,
    kappa::Float64,
    cp::Float64,
    lv::Float64,
    rv::Float64,
)
    pressure_b <= pressure_a ||
        throw(ArgumentError("moist parcel ascent requires nonincreasing pressure"))
    delta_log_pressure = log(pressure_b / pressure_a)
    derivative_a = _bm_moist_derivative(temperature_a, mixing_ratio_a, kappa, cp, lv, rv)
    temperature_mid = temperature_a + 0.5 * derivative_a * delta_log_pressure
    pressure_mid = 0.5 * (pressure_a + pressure_b)
    mixing_ratio_mid = Saturation_Mixing_Ratio(temperature_mid, pressure_mid, epsilon)
    derivative_mid =
        _bm_moist_derivative(temperature_mid, mixing_ratio_mid, kappa, cp, lv, rv)
    temperature_b = temperature_a + derivative_mid * delta_log_pressure
    mixing_ratio_b = Saturation_Mixing_Ratio(temperature_b, pressure_b, epsilon)
    return temperature_b, mixing_ratio_b
end

function _betts_miller_column!(
    work::Betts_Miller_Work,
    state::Betts_Miller_State,
    temperature::AbstractVector{<:Real},
    humidity::AbstractVector{<:Real},
    p_full::AbstractVector{<:Real},
    p_half::AbstractVector{<:Real},
    rd::Float64,
    rv::Float64,
    cp::Float64,
    lv::Float64,
    grav::Float64,
    kappa::Float64,
    use_virtual_temperature::Bool,
)
    _validate_bm_column(temperature, humidity, p_full, p_half)
    nd = state.nd
    length(p_full) == nd || throw(DimensionMismatch("column does not match BM state"))

    tp = work.parcel_temperature
    rp = work.parcel_mixing_ratio
    rs = work.parcel_saturation_mixing_ratio
    tref = work.reference_temperature
    qref = work.reference_humidity
    tdot = work.temperature_tendency
    qdot = work.humidity_tendency

    # Treat roundoff-scale negative inputs as dry without mutating the caller.
    tp .= temperature
    @. qref = max(humidity, 0.0)
    @. rp = qref / (1.0 - qref)
    fill!(rs, 0.0)
    tref .= temperature
    fill!(tdot, 0.0)
    fill!(qdot, 0.0)

    epsilon = rd / rv
    virtual_coefficient = use_virtual_temperature ? rv / rd - 1.0 : 0.0
    pstar = 1.0e5
    surface = nd
    t0 = Float64(temperature[surface])
    r0 = rp[surface]
    rs0 = Saturation_Mixing_Ratio(t0, Float64(p_full[surface]), epsilon)

    cape = 0.0
    cin = 0.0
    lcl = 0
    lfc = 0
    lzb = 0
    first_buoyant = true

    if r0 >= rs0
        lcl = surface
        tp[surface] = t0 + (r0 - rs0) / (cp / lv + lv * rs0 / (rv * t0^2))
        rp[surface] =
            Saturation_Mixing_Ratio(tp[surface], Float64(p_full[surface]), epsilon)
        rs[surface] = rp[surface]
    else
        theta0 = t0 * (pstar / Float64(p_full[surface]))^kappa
        condensation = _bm_lcl(
            theta0,
            r0,
            Float64(p_full[1]),
            Float64(p_full[surface]),
            epsilon,
            kappa,
        )
        if !condensation.found
            # The parcel never condenses in the represented domain. Retain its
            # dry ascent diagnostics without inventing a moist adjustment.
            for k = surface:-1:1
                tp[k] = t0 * (Float64(p_full[k]) / Float64(p_full[surface]))^kappa
                rp[k] = r0
                rs[k] = Saturation_Mixing_Ratio(tp[k], Float64(p_full[k]), epsilon)
                buoyancy = _bm_buoyancy(
                    tp[k], rp[k], temperature[k], humidity[k], virtual_coefficient,
                )
                cin +=
                    rd * max(-buoyancy, 0.0) *
                    _bm_layer_log_pressure(p_full, p_half, k)
            end
            return (;
                active = false,
                regime = :none,
                adjustment_top = 0,
                top_fraction = 0.0,
                lcl = 0,
                lfc = 0,
                lzb = 0,
                cape = 0.0,
                cin,
                precipitation = 0.0,
            )
        end
        plcl, tlcl = condensation.pressure, condensation.temperature

        k = surface
        while k >= 1 && p_full[k] > plcl && !isapprox(p_full[k], plcl; rtol = 32 * eps(Float64))
            tp[k] = t0 * (Float64(p_full[k]) / Float64(p_full[surface]))^kappa
            rp[k] = r0
            rs[k] = Saturation_Mixing_Ratio(tp[k], Float64(p_full[k]), epsilon)
            buoyancy = _bm_buoyancy(
                tp[k], rp[k], temperature[k], humidity[k], virtual_coefficient,
            )
            cin -= rd * buoyancy * _bm_layer_log_pressure(p_full, p_half, k)
            k -= 1
        end
        lcl = k
        if isapprox(p_full[lcl], plcl; rtol = 32 * eps(Float64))
            # Snap a condensation level coincident with a full level to that
            # level, so bisection roundoff cannot move it into an adjacent cell.
            tp[lcl] = t0 * (Float64(p_full[lcl]) / Float64(p_full[surface]))^kappa
            rp[lcl] = Saturation_Mixing_Ratio(tp[lcl], Float64(p_full[lcl]), epsilon)
        else
            tp[lcl], rp[lcl] =
                _bm_moist_rk2(tlcl, r0, plcl, Float64(p_full[lcl]), epsilon, kappa, cp, lv, rv)
        end
        rs[lcl] = rp[lcl]
        if tp[lcl] < BM_MIN_PARCEL_TEMPERATURE
            return (;
                active = false,
                regime = :none,
                adjustment_top = 0,
                top_fraction = 0.0,
                lcl,
                lfc = 0,
                lzb = 0,
                cape = 0.0,
                cin = 0.0,
                precipitation = 0.0,
            )
        end

        layer_factor = _bm_layer_log_pressure(p_full, p_half, lcl)
        buoyancy = _bm_buoyancy(
            tp[lcl], rp[lcl], temperature[lcl], humidity[lcl], virtual_coefficient,
        )
        if buoyancy < 0
            cin -= rd * buoyancy * layer_factor
        elseif buoyancy > 0
            cape += rd * buoyancy * layer_factor
            first_buoyant = false
            lfc = lcl
        end
    end

    for k = lcl-1:-1:1
        tp[k], rp[k] = _bm_moist_rk2(
            tp[k+1],
            rp[k+1],
            Float64(p_full[k+1]),
            Float64(p_full[k]),
            epsilon,
            kappa,
            cp,
            lv,
            rv,
        )
        rs[k] = rp[k]
        if tp[k] < BM_MIN_PARCEL_TEMPERATURE && first_buoyant
            return (;
                active = false,
                regime = :none,
                adjustment_top = 0,
                top_fraction = 0.0,
                lcl,
                lfc = 0,
                lzb = 0,
                cape = 0.0,
                cin = 0.0,
                precipitation = 0.0,
            )
        end

        layer_factor = _bm_layer_log_pressure(p_full, p_half, k)
        buoyancy = _bm_buoyancy(tp[k], rp[k], temperature[k], humidity[k], virtual_coefficient)
        if buoyancy < 0
            if first_buoyant
                cin -= rd * buoyancy * layer_factor
            else
                lzb = k + 1
                break
            end
        elseif buoyancy > 0
            cape += rd * buoyancy * layer_factor
            if first_buoyant
                first_buoyant = false
                lfc = k
            end
        end
    end

    if first_buoyant || cape <= 0
        fill!(tdot, 0.0)
        fill!(qdot, 0.0)
        return (;
            active = false,
            regime = :none,
            adjustment_top = 0,
            top_fraction = 0.0,
            lcl,
            lfc = 0,
            lzb = 0,
            cape = 0.0,
            cin = 0.0,
            precipitation = 0.0,
        )
    end
    lzb == 0 && (lzb = 1)

    for k = lzb:surface
        tref[k] = tp[k]
        # Relative humidity scales vapor pressure, not mixing ratio. Recover
        # saturation vapor pressure from the separately stored saturation ratio.
        qref[k] = _bm_reference_humidity(rs[k], Float64(p_full[k]), state.relative_humidity, epsilon)
        tdot[k] = (tref[k] - Float64(temperature[k])) / state.tau
        qdot[k] = (qref[k] - max(Float64(humidity[k]), 0.0)) / state.tau
    end

    adjustment = _bm_close_adjustment!(
        work, state, temperature, humidity, p_half, lzb, surface, cp, lv, grav,
    )
    return merge(adjustment, (; lcl, lfc, lzb, cape, cin))
end

"""
Run the Betts-Miller calculation for one top-to-bottom atmospheric column.

`lcl = 0` indicates no condensation in the represented ascent domain.
`parcel_mixing_ratio` is the actual parcel mixing ratio: it is conserved below
the LCL and saturated above it. `parcel_saturation_mixing_ratio` stores the
saturation mixing ratio used to construct reference humidity at visited levels.
Levels above the end of ascent retain the input temperature and mixing ratio.
`regime` is `:none`, `:deep`, or `:shallow`; `active` includes zero-rain shallow
transport. `adjustment_top` and `top_fraction` describe the accepted adjustment,
while `lzb` retains the diagnosed buoyancy limit. Inactive columns use top 0
and fraction 0. Shallow references are effective full-cell targets after
fractional penetration and the temperature correction.
"""
function Betts_Miller_Column(
    state::Betts_Miller_State,
    atmo_data::Atmo_Data,
    temperature::AbstractVector{<:Real},
    humidity::AbstractVector{<:Real},
    p_full::AbstractVector{<:Real},
    p_half::AbstractVector{<:Real},
)
    work = Betts_Miller_Work(state.nd)
    diagnostics = _betts_miller_column!(
        work,
        state,
        temperature,
        humidity,
        p_full,
        p_half,
        atmo_data.rdgas,
        atmo_data.rvgas,
        atmo_data.cp_air,
        atmo_data.Lv,
        atmo_data.grav,
        atmo_data.kappa,
        atmo_data.use_virtual_temperature,
    )
    return merge(
        diagnostics,
        (
            temperature_tendency = copy(work.temperature_tendency),
            humidity_tendency = copy(work.humidity_tendency),
            parcel_temperature = copy(work.parcel_temperature),
            parcel_mixing_ratio = copy(work.parcel_mixing_ratio),
            parcel_saturation_mixing_ratio = copy(work.parcel_saturation_mixing_ratio),
            reference_temperature = copy(work.reference_temperature),
            reference_humidity = copy(work.reference_humidity),
        ),
    )
end

"""
    Betts_Miller!(state, atmo_data, T, q, p_full, p_half, bm_dt, bm_dq, bm_precip)

Calculate Betts-Miller temperature and humidity tendencies and precipitation
rate on every grid column. Output buffers are overwritten.
"""
function Betts_Miller!(
    state::Betts_Miller_State,
    atmo_data::Atmo_Data,
    temperature::Array{Float64,3},
    humidity::Array{Float64,3},
    p_full::Array{Float64,3},
    p_half::Array{Float64,3},
    bm_temperature_tendency::Array{Float64,3},
    bm_humidity_tendency::Array{Float64,3},
    bm_precipitation::Array{Float64,3},
)
    size(temperature) == size(humidity) == size(p_full) ||
        throw(DimensionMismatch("Betts-Miller full-level fields must have equal sizes"))
    size(temperature) == size(bm_temperature_tendency) == size(bm_humidity_tendency) ||
        throw(DimensionMismatch("Betts-Miller tendency fields have incorrect sizes"))
    nλ, nθ, nd = size(temperature)
    size(p_half) == (nλ, nθ, nd + 1) ||
        throw(DimensionMismatch("Betts-Miller p_half has incorrect size"))
    size(bm_precipitation) == (nλ, nθ, 1) ||
        throw(DimensionMismatch("Betts-Miller precipitation has incorrect size"))
    nd == state.nd || throw(DimensionMismatch("Betts-Miller state has incorrect nd"))

    fill!(bm_temperature_tendency, 0.0)
    fill!(bm_humidity_tendency, 0.0)
    fill!(bm_precipitation, 0.0)

    @threads for j = 1:nθ
        work = state.work[Threads.threadid()]
        for i = 1:nλ
            diagnostics = _betts_miller_column!(
                work,
                state,
                @view(temperature[i, j, :]),
                @view(humidity[i, j, :]),
                @view(p_full[i, j, :]),
                @view(p_half[i, j, :]),
                atmo_data.rdgas,
                atmo_data.rvgas,
                atmo_data.cp_air,
                atmo_data.Lv,
                atmo_data.grav,
                atmo_data.kappa,
                atmo_data.use_virtual_temperature,
            )
            @views bm_temperature_tendency[i, j, :] .= work.temperature_tendency
            @views bm_humidity_tendency[i, j, :] .= work.humidity_tendency
            bm_precipitation[i, j, 1] = diagnostics.precipitation
        end
    end
    return nothing
end
