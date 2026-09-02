module Initial_Conditions

using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Dyn_Data_Module
using ..Experiment_Configuration
using ..Variable_Mappings_Module
using ..Vert_Coordinate_Module
using ..Press_And_Geopot_Module: Pressure_Variables!
using ..Spectral_Dynamics_Module: Get_Topography!
using ..Restart_Manager_Module

export Initialize_Atmos_State!



"""
    Initialize_Atmos_State!(mesh, atmo_data, dyn_data, ic_source)

Acts as a bridge between the user's initial conditions and the model's internal state.
- If `ic_source` is a **Function**: Calls `ic_source(mesh, atmo_data, dyn_data)`.
- If `ic_source` is a **Dict/NamedTuple**: Automatically transforms grid arrays to spectral coefficients.
- If `ic_source` is a **Symbol**: Dispatches to built-in benchmark cases.
"""
function Initialize_Atmos_State!(
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    vert_coord::Union{Vert_Coordinate,Nothing},
    config::Model_Config,
)

    # --- Restart --- #
    if config.is_restart
        if isfile(config.restart_file)
            Load_Restart_File!(dyn_data, config.restart_file)
            @info "Initialization complete: Loaded warm start from $(config.restart_file)"
            return # <--- Crucial: Stop here!
        else
            error("Restart requested but file not found: $(config.restart_file)")
        end
    end
    # --- Restart --- #

    # --- Cold start --- #
    ic_source = config.initial_condition

    if ic_source isa Function
        # 1. Custom Function provided by user
        ic_source(mesh, atmo_data, dyn_data, vert_coord)

    elseif ic_source isa Union{Dict,NamedTuple}
        # 2. Data Injection (Arrays)
        # We pass model_type to handle specific aliases (like h -> ps in SW)
        Load_From_Arrays!(mesh, dyn_data, ic_source, config.model_type)

    elseif ic_source isa Symbol
        # 3. Standard Benchmarks
        Dispatch_Standard_IC!(mesh, atmo_data, dyn_data, vert_coord, ic_source, config)

    else
        error("Unsupported Initial Condition Type: $(typeof(ic_source))")
    end
    # --- Cold start --- #

end



"""
    Load_From_Arrays!(mesh, dyn_data, data)

Helper that takes user-provided grid arrays (u, v, t, ps, etc.) 
and handles the spectral transformations automatically.
"""
function Load_From_Arrays!(mesh, dyn_data, data, model_type::Symbol)

    # Get the authoritative map for this model type
    var_map = Get_Dyn_Var_Map(dyn_data, Val(model_type))

    # Iterate and Load
    for (key, user_array) in data
        if key == :u || key == :v
            continue
        end # Handled later

        if haskey(var_map, key)
            target_field = var_map[key]

            if target_field isa AbstractArray
                if size(target_field) == size(user_array)
                    target_field .= user_array
                else
                    error(
                        "Dimension mismatch for :$key. Expected $(size(target_field)), got $(size(user_array))",
                    )
                end

                # Transform Prognostics (Grid -> Spectral)
                if key == :t
                    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_t_c, dyn_data.spe_t_c)
                elseif key == :lnps || key == :ps || key == :h
                    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_lnps, dyn_data.spe_lnps_c)
                end
            end
        end
    end

    # Handle Winds Explicitly
    if haskey(data, :u) && haskey(data, :v)
        dyn_data.grid_u_c .= data[:u]
        dyn_data.grid_v_c .= data[:v]

        Vor_Div_From_Grid_UV!(
            mesh,
            dyn_data.grid_u_c,
            dyn_data.grid_v_c,
            dyn_data.spe_vor_c,
            dyn_data.spe_div_c,
        )

        Trans_Spherical_To_Grid!(mesh, dyn_data.spe_vor_c, dyn_data.grid_vor)
        Trans_Spherical_To_Grid!(mesh, dyn_data.spe_div_c, dyn_data.grid_div)
    end
end



"""
    Dispatch_Standard_IC!

Legacy support for named benchmarks.
"""
function Dispatch_Standard_IC!(
    mesh,
    atmo_data,
    dyn_data,
    vert_coord,
    ic_name::Symbol,
    config,
)

    if ic_name == :Barotropic_Jet
        Init_Barotropic_Jet!(mesh, dyn_data)

    elseif ic_name == :Shallow_Water_Test
        Init_Shallow_Water_Test!(mesh, dyn_data, config)

    elseif ic_name == :Moist_Spinup
        Init_Moist_Aquaplanet!(mesh, atmo_data, dyn_data, vert_coord, config)

    else
        error("Unknown Benchmark IC: $ic_name")
    end
end

# ==============================================================================
# 1. Barotropic Initialization
#    (Logic extracted from exp/Barotropic/Barotropic.jl)
# ==============================================================================
function Init_Barotropic_Jet!(mesh::Spectral_Spherical_Mesh, dyn_data::Dyn_Data)

    grid_u = dyn_data.grid_u_c
    grid_v = dyn_data.grid_v_c
    cosθ, sinθ = mesh.cosθ, mesh.sinθ

    # 1. Set Zonal Wind Profile (Winter-like Jet)
    # u = 25cosθ - 30cosθ³ + 300sin²θcos⁶θ
    for k = 1:mesh.nd
        for j = 1:mesh.nθ
            for i = 1:mesh.nλ
                c_theta = cosθ[j]
                s_theta = sinθ[j]
                grid_u[i, j, k] =
                    25.0 * c_theta - 30.0 * c_theta^3 + 300.0 * s_theta^2 * c_theta^6
                grid_v[i, j, k] = 0.0
            end
        end
    end

    # 2. Convert to Spectral Space to get smooth Vor/Div
    spe_vor_c, spe_div_c = dyn_data.spe_vor_c, dyn_data.spe_div_c
    Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)
    Trans_Spherical_To_Grid!(mesh, spe_vor_c, dyn_data.grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_c, dyn_data.grid_div)

    # 3. Add Gaussian Perturbation to Vorticity
    # Params: m=4, lat=45deg, width=15deg
    m, θ0, θw, A = 4.0, 45.0 * pi / 180, 15.0 * pi / 180.0, 8.0e-5

    grid_vor = dyn_data.grid_vor
    λc, θc = mesh.λc, mesh.θc

    for k = 1:mesh.nd
        for j = 1:mesh.nθ
            for i = 1:mesh.nλ
                # A/2 * cos(lat) * exp(...) * cos(m*lon)
                perturbation =
                    (A / 2.0) * cosθ[j] * exp(-((θc[j] - θ0) / θw)^2) * cos(m * λc[i])
                grid_vor[i, j, k] += perturbation
            end
        end
    end

    # 4. Update Spectral Coefficients with Perturbation
    Trans_Grid_To_Spherical!(mesh, grid_vor, spe_vor_c)

    # Optional: Re-compute U/V from perturbed vorticity if strict consistency is needed,
    # but for barotropic evolution, Vor is the prognostic variable.
end

# ==============================================================================
# 2. Shallow Water Initialization
#    (Logic extracted from exp/Shallow_Water/SW.jl)
# ==============================================================================
function Init_Shallow_Water_Test!(
    mesh::Spectral_Spherical_Mesh,
    dyn_data::Dyn_Data,
    config::Model_Config,
)

    params = config.physics_params
    h_0 = get(params, "h_0", 3.0e4)
    h_amp = get(params, "h_amp", 2.0e4)
    h_lon = get(params, "h_lon", 90.0) * pi / 180
    h_lat = get(params, "h_lat", 25.0) * pi / 180
    h_width = get(params, "h_width", 15.0) * pi / 180
    h_itcz = get(params, "h_itcz", 1.0e5)
    itcz_width = get(params, "itcz_width", 4.0) * pi / 180

    # 1. Initialize U, V to zero
    dyn_data.grid_u_c .= 0.0
    dyn_data.grid_v_c .= 0.0
    Vor_Div_From_Grid_UV!(
        mesh,
        dyn_data.grid_u_c,
        dyn_data.grid_v_c,
        dyn_data.spe_vor_c,
        dyn_data.spe_div_c,
    )
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_vor_c, dyn_data.grid_vor)
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_div_c, dyn_data.grid_div)

    # 2. Initialize Height Field (Prognostic Variable)
    dyn_data.grid_lnps .= h_0
    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_lnps, dyn_data.spe_lnps_c)

    # 3. Construct Equilibrium Height (Forcing Field)
    # Stored in grid_geopots for the physics routine
    h_eq = dyn_data.grid_geopots
    λc, θc = mesh.λc, mesh.θc

    for k = 1:mesh.nd
        for j = 1:mesh.nθ
            d2 = (θc[j] / itcz_width)^2
            for i = 1:mesh.nλ
                x2 = ((λc[i] - h_lon) / (2.0 * h_width))^2
                y2 = ((θc[j] - h_lat) / h_width)^2

                h_eq[i, j, k] = h_0 + h_amp * exp(-(x2 + y2)) + h_itcz * exp(-d2)
            end
        end
    end
end

# ==============================================================================
# 3. Moist Held--Suarez aquaplanet initialization
#
# Thatcher and Jablonowski (2016), Appendix A, use the shallow-atmosphere
# balanced state of Ullrich et al. (2014), with their analytic humidity profile.
# Pressure-coordinate models must first invert the balanced pressure relation to
# obtain the height corresponding to each (latitude, pressure) point.
# ==============================================================================

const MOIST_AQUAPLANET_P0 = 1.0e5
const MOIST_AQUAPLANET_GAMMA = 0.005
const MOIST_AQUAPLANET_B = 2.0
const MOIST_AQUAPLANET_K = 3.0
const MOIST_AQUAPLANET_T_EQUATOR = 310.0
const MOIST_AQUAPLANET_T_POLE = 240.0

"""Return the Ullrich shallow-atmosphere basic state at `(latitude, pressure)`."""
function Ullrich_Shallow_Basic_State(
    atmo_data::Atmo_Data,
    latitude::Float64,
    pressure::Float64,
)
    pressure > 0.0 || throw(DomainError(pressure, "initial pressure must be positive"))

    rdgas = atmo_data.rdgas
    grav = atmo_data.grav
    radius = atmo_data.radius
    omega = atmo_data.omega

    Γ = MOIST_AQUAPLANET_GAMMA
    b = MOIST_AQUAPLANET_B
    k = MOIST_AQUAPLANET_K
    te = MOIST_AQUAPLANET_T_EQUATOR
    tp = MOIST_AQUAPLANET_T_POLE
    t0 = 0.5 * (te + tp)
    scale_height = rdgas * t0 / grav

    A = 1.0 / Γ
    B = (te - tp) / ((te + tp) * tp)
    C = 0.5 * (k + 2) * (te - tp) / (te * tp)

    coslat = cos(latitude)
    latitude_factor = coslat^k - (k / (k + 2)) * coslat^(k + 2)

    # Ullrich et al. (2014), Appendix C.2, equations (40), (43), and
    # (44). A hypsometric estimate is a more reliable initial guess than the
    # fixed 10 km value suggested in the paper.
    z = max(0.0, -(rdgas * t0 / grav) * log(pressure / MOIST_AQUAPLANET_P0))
    for _ = 1:12
        scaled_z = z / (b * scale_height)
        gaussian = exp(-scaled_z^2)
        integral_tau1 = A * (exp(Γ * z / t0) - 1.0) + B * z * gaussian
        integral_tau2 = C * z * gaussian
        residual =
            log(pressure / MOIST_AQUAPLANET_P0) +
            (grav / rdgas) * (integral_tau1 - integral_tau2 * latitude_factor)

        tau1 =
            (1.0 / t0) * exp(Γ * z / t0) +
            B * (1.0 - 2.0 * scaled_z^2) * gaussian
        tau2 = C * (1.0 - 2.0 * scaled_z^2) * gaussian
        derivative = (grav / rdgas) * (tau1 - tau2 * latitude_factor)
        step = residual / derivative
        z -= step
        abs(step) <= 1.0e-9 * max(1.0, z) && break
    end

    scaled_z = z / (b * scale_height)
    gaussian = exp(-scaled_z^2)
    tau1 =
        (1.0 / t0) * exp(Γ * z / t0) +
        B * (1.0 - 2.0 * scaled_z^2) * gaussian
    tau2 = C * (1.0 - 2.0 * scaled_z^2) * gaussian
    virtual_temperature = 1.0 / (tau1 - tau2 * latitude_factor)

    integral_tau2 = C * z * gaussian
    wind_proxy =
        (grav * k / radius) *
        integral_tau2 *
        (coslat^(k - 1) - coslat^(k + 1)) *
        virtual_temperature
    rotation_speed = omega * radius * coslat
    zonal_wind =
        -rotation_speed +
        sqrt(max(0.0, rotation_speed^2 + radius * coslat * wind_proxy))

    return (; height=z, virtual_temperature, zonal_wind)
end

"""Return the rotational wind perturbation from Ullrich et al. (2014)."""
function Ullrich_Wind_Perturbation(
    atmo_data::Atmo_Data,
    longitude::Float64,
    latitude::Float64,
    height::Float64,
)
    perturbation_top = 1.5e4
    (0.0 <= height <= perturbation_top) || return (u=0.0, v=0.0)

    center_lon = pi / 9.0
    center_lat = 2.0 * pi / 9.0
    angular_radius = 1.0 / 6.0
    amplitude = 1.0

    central_angle = acos(
        clamp(
            sin(center_lat) * sin(latitude) +
            cos(center_lat) * cos(latitude) * cos(longitude - center_lon),
            -1.0,
            1.0,
        ),
    )
    central_angle <= angular_radius || return (u=0.0, v=0.0)
    abs(sin(central_angle)) > 1.0e-14 || return (u=0.0, v=0.0)

    vertical_fraction = height / perturbation_top
    taper = 1.0 - 3.0 * vertical_fraction^2 + 2.0 * vertical_fraction^3
    horizontal_phase = pi * central_angle / (2.0 * angular_radius)
    factor =
        (16.0 * amplitude / (3.0 * sqrt(3.0))) *
        taper *
        cos(horizontal_phase)^3 *
        sin(horizontal_phase) /
        sin(central_angle)

    u =
        -factor *
        (-sin(center_lat) * cos(latitude) +
         cos(center_lat) * sin(latitude) * cos(longitude - center_lon))
    v = factor * cos(center_lat) * sin(longitude - center_lon)
    return (; u, v)
end

"""
    Init_Moist_Aquaplanet!(mesh, atmo_data, dyn_data, vert_coord, config)

Initialize the moist Held--Suarez experiment from the balanced shallow-atmosphere
baroclinic-wave state prescribed in Appendix A of Thatcher and Jablonowski
(2016). The dry analytic temperature is virtual temperature; actual temperature
is therefore obtained using the model's own water-vapor gas constant.
"""
function Init_Moist_Aquaplanet!(
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    vert_coord::Vert_Coordinate,
    config::Model_Config,
)
    Get_Topography!(dyn_data.grid_geopots)

    dyn_data.grid_lnps .= log(MOIST_AQUAPLANET_P0)
    dyn_data.grid_ps_c .= MOIST_AQUAPLANET_P0
    Pressure_Variables!(
        vert_coord,
        dyn_data.grid_ps_c,
        dyn_data.grid_p_half,
        dyn_data.grid_Δp,
        dyn_data.grid_lnp_half,
        dyn_data.grid_p_full,
        dyn_data.grid_lnp_full,
    )

    moisture_coefficient = atmo_data.rvgas / atmo_data.rdgas - 1.0
    humidity_floor = Float64(get(config.physics_params, "initial_humidity_floor", 0.0))
    0.0 <= humidity_floor < 1.0 ||
        throw(ArgumentError("initial_humidity_floor must satisfy 0 <= q < 1"))

    for k = 1:mesh.nd, j = 1:mesh.nθ, i = 1:mesh.nλ
        longitude = mesh.λc[i]
        latitude = mesh.θc[j]
        pressure = dyn_data.grid_p_full[i, j, k]
        basic_state = Ullrich_Shallow_Basic_State(atmo_data, latitude, pressure)
        perturbation =
            Ullrich_Wind_Perturbation(atmo_data, longitude, latitude, basic_state.height)

        q = 0.0
        if config.moisture_processes
            if pressure >= 1.0e4
                q =
                    0.018 *
                    exp(-(latitude / (2.0 * pi / 9.0))^4) *
                    exp(-(((pressure / MOIST_AQUAPLANET_P0 - 1.0) *
                           MOIST_AQUAPLANET_P0 / 3.0e4)^2))
                q = max(q, humidity_floor)
            end
        end

        dyn_data.grid_q_c[i, j, k] = q
        dyn_data.grid_u_c[i, j, k] = basic_state.zonal_wind + perturbation.u
        dyn_data.grid_v_c[i, j, k] = perturbation.v
    end

    # The dry dynamical fields remain spectrally truncated. Specific humidity is
    # grid-only and retains the analytic Thatcher--Jablonowski initial profile.
    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_lnps, dyn_data.spe_lnps_c)
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_lnps_c, dyn_data.grid_lnps)
    dyn_data.grid_ps_c .= exp.(dyn_data.grid_lnps)

    if !config.moisture_processes
        dyn_data.grid_q_c .= 0.0
    end

    # Convert the prescribed virtual temperature with the same grid humidity
    # used by the moist equation of state.
    for k = 1:mesh.nd, j = 1:mesh.nθ, i = 1:mesh.nλ
        basic_state = Ullrich_Shallow_Basic_State(
            atmo_data,
            mesh.θc[j],
            dyn_data.grid_p_full[i, j, k],
        )
        dyn_data.grid_t_c[i, j, k] =
            basic_state.virtual_temperature /
            (1.0 + moisture_coefficient * dyn_data.grid_q_c[i, j, k])
    end
    Trans_Grid_To_Spherical!(mesh, dyn_data.grid_t_c, dyn_data.spe_t_c)
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_t_c, dyn_data.grid_t_c)

    Vor_Div_From_Grid_UV!(
        mesh,
        dyn_data.grid_u_c,
        dyn_data.grid_v_c,
        dyn_data.spe_vor_c,
        dyn_data.spe_div_c,
    )
    UV_Grid_From_Vor_Div!(
        mesh,
        dyn_data.spe_vor_c,
        dyn_data.spe_div_c,
        dyn_data.grid_u_c,
        dyn_data.grid_v_c,
    )
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_vor_c, dyn_data.grid_vor)
    Trans_Spherical_To_Grid!(mesh, dyn_data.spe_div_c, dyn_data.grid_div)

    dyn_data.spe_vor_p .= dyn_data.spe_vor_c
    dyn_data.spe_div_p .= dyn_data.spe_div_c
    dyn_data.spe_lnps_p .= dyn_data.spe_lnps_c
    dyn_data.spe_t_p .= dyn_data.spe_t_c
    dyn_data.grid_u_p .= dyn_data.grid_u_c
    dyn_data.grid_v_p .= dyn_data.grid_v_c
    dyn_data.grid_ps_p .= dyn_data.grid_ps_c
    dyn_data.grid_t_p .= dyn_data.grid_t_c
    dyn_data.grid_q_p .= dyn_data.grid_q_c
end

end
