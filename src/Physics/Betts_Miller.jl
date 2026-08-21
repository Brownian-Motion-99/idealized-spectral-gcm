using Base.Threads

mutable struct Betts_Miller_Work
    parcel_temperature::Vector{Float64}
    parcel_mixing_ratio::Vector{Float64}
    reference_temperature::Vector{Float64}
    reference_humidity::Vector{Float64}
    temperature_tendency::Vector{Float64}
    humidity_tendency::Vector{Float64}
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
    )
end

"""
    Betts_Miller_State(nd; tau=7200.0, relative_humidity=0.8)

Time-independent configuration and per-thread work storage for the Betts-Miller
convective adjustment. The scheme returns rates; the model time step is not part
of the column calculation.
"""
struct Betts_Miller_State
    tau::Float64
    relative_humidity::Float64
    nd::Int
    work::Vector{Betts_Miller_Work}
end

function Betts_Miller_State(nd::Int; tau::Real = 7200.0, relative_humidity::Real = 0.8)
    tau = Float64(tau)
    relative_humidity = Float64(relative_humidity)
    isfinite(tau) && tau > 0 || throw(ArgumentError("bm_tau must be positive and finite"))
    isfinite(relative_humidity) && 0 < relative_humidity <= 1 ||
        throw(ArgumentError("bm_relative_humidity must lie in (0, 1]"))
    return Betts_Miller_State(
        tau,
        relative_humidity,
        nd,
        [Betts_Miller_Work(nd) for _ = 1:Threads.nthreads()],
    )
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
        _bm_saturation_mixing_ratio(parcel_temperature, pressure, epsilon) - r0
    end

    lo = log(p_top)
    hi = log(p_surface)
    if residual(lo) >= 0
        pressure = p_top
        return pressure, theta0 * (pressure / pstar)^kappa
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
    return pressure, theta0 * (pressure / pstar)^kappa
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
    delta_log_pressure = log(pressure_b / pressure_a)
    derivative_a = _bm_moist_derivative(temperature_a, mixing_ratio_a, kappa, cp, lv, rv)
    temperature_mid = temperature_a + 0.5 * derivative_a * delta_log_pressure
    pressure_mid = 0.5 * (pressure_a + pressure_b)
    mixing_ratio_mid = _bm_saturation_mixing_ratio(temperature_mid, pressure_mid, epsilon)
    derivative_mid =
        _bm_moist_derivative(temperature_mid, mixing_ratio_mid, kappa, cp, lv, rv)
    temperature_b = temperature_a + derivative_mid * delta_log_pressure
    mixing_ratio_b = _bm_saturation_mixing_ratio(temperature_b, pressure_b, epsilon)
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
)
    _validate_bm_column(temperature, humidity, p_full, p_half)
    nd = state.nd
    length(p_full) == nd || throw(DimensionMismatch("column does not match BM state"))

    tp = work.parcel_temperature
    rp = work.parcel_mixing_ratio
    tref = work.reference_temperature
    qref = work.reference_humidity
    tdot = work.temperature_tendency
    qdot = work.humidity_tendency

    # Spherical-harmonic transforms can create tiny negative grid-point
    # undershoots even when the spectral humidity field is physically valid.
    # Diagnose those points as dry without mutating the prognostic state.
    tp .= temperature
    @. qref = max(humidity, 0.0)
    @. rp = qref / (1.0 - qref)
    tref .= temperature
    fill!(tdot, 0.0)
    fill!(qdot, 0.0)

    epsilon = rd / rv
    pstar = 1.0e5
    surface = nd
    t0 = Float64(temperature[surface])
    r0 = rp[surface]
    rs0 = _bm_saturation_mixing_ratio(t0, Float64(p_full[surface]), epsilon)

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
            _bm_saturation_mixing_ratio(tp[surface], Float64(p_full[surface]), epsilon)
    elseif r0 > 0
        theta0 = t0 * (pstar / Float64(p_full[surface]))^kappa
        plcl, tlcl = _bm_lcl(
            theta0,
            r0,
            Float64(p_full[1]),
            Float64(p_full[surface]),
            epsilon,
            kappa,
        )

        k = surface
        while k >= 1 && p_full[k] > plcl
            tp[k] = theta0 * (Float64(p_full[k]) / pstar)^kappa
            rp[k] = _bm_saturation_mixing_ratio(tp[k], Float64(p_full[k]), epsilon)
            cin +=
                rd *
                (Float64(temperature[k]) - tp[k]) *
                _bm_layer_log_pressure(p_full, p_half, k)
            k -= 1
        end
        lcl = max(k, 2)
        tp[lcl], rp[lcl] =
            _bm_moist_rk2(tlcl, r0, plcl, Float64(p_full[lcl]), epsilon, kappa, cp, lv, rv)
        if tp[lcl] < BM_MIN_PARCEL_TEMPERATURE
            return (;
                active = false,
                lcl,
                lfc = 0,
                lzb = 0,
                cape = 0.0,
                cin = 0.0,
                precipitation = 0.0,
            )
        end

        layer_factor = _bm_layer_log_pressure(p_full, p_half, lcl)
        if tp[lcl] < temperature[lcl]
            cin += rd * (Float64(temperature[lcl]) - tp[lcl]) * layer_factor
        else
            cape += rd * (tp[lcl] - Float64(temperature[lcl])) * layer_factor
            first_buoyant = false
            lfc = lcl
        end
    else
        return (;
            active = false,
            lcl = 0,
            lfc = 0,
            lzb = 0,
            cape = 0.0,
            cin = 0.0,
            precipitation = 0.0,
        )
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
        if tp[k] < BM_MIN_PARCEL_TEMPERATURE && first_buoyant
            return (;
                active = false,
                lcl,
                lfc = 0,
                lzb = 0,
                cape = 0.0,
                cin = 0.0,
                precipitation = 0.0,
            )
        end

        layer_factor = _bm_layer_log_pressure(p_full, p_half, k)
        if tp[k] < temperature[k]
            if first_buoyant
                cin += rd * (Float64(temperature[k]) - tp[k]) * layer_factor
            else
                lzb = k + 1
                break
            end
        else
            cape += rd * (tp[k] - Float64(temperature[k])) * layer_factor
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
        reference_mixing_ratio = state.relative_humidity * rp[k]
        qref[k] = reference_mixing_ratio / (1.0 + reference_mixing_ratio)
        tdot[k] = (tref[k] - Float64(temperature[k])) / state.tau
        qdot[k] = (qref[k] - max(Float64(humidity[k]), 0.0)) / state.tau
    end

    moisture_precipitation = 0.0
    thermal_precipitation = 0.0
    for k = lzb:surface
        layer_mass = (Float64(p_half[k+1]) - Float64(p_half[k])) / grav
        moisture_precipitation -= qdot[k] * layer_mass
        thermal_precipitation += (cp / lv) * tdot[k] * layer_mass
    end

    if moisture_precipitation <= 0 || thermal_precipitation <= 0
        fill!(tdot, 0.0)
        fill!(qdot, 0.0)
        return (; active = false, lcl, lfc, lzb, cape, cin, precipitation = 0.0)
    elseif moisture_precipitation > thermal_precipitation
        scale = thermal_precipitation / moisture_precipitation
        @views qdot[lzb:surface] .*= scale
        precipitation = thermal_precipitation
    else
        scale = moisture_precipitation / thermal_precipitation
        @views tdot[lzb:surface] .*= scale
        precipitation = moisture_precipitation
    end

    return (; active = true, lcl, lfc, lzb, cape, cin, precipitation)
end

"""Run the Betts-Miller calculation for one top-to-bottom atmospheric column."""
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
    )
    return merge(
        diagnostics,
        (
            temperature_tendency = copy(work.temperature_tendency),
            humidity_tendency = copy(work.humidity_tendency),
            parcel_temperature = copy(work.parcel_temperature),
            parcel_mixing_ratio = copy(work.parcel_mixing_ratio),
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
            )
            @views bm_temperature_tendency[i, j, :] .= work.temperature_tendency
            @views bm_humidity_tendency[i, j, :] .= work.humidity_tendency
            bm_precipitation[i, j, 1] = diagnostics.precipitation
        end
    end
    return nothing
end
