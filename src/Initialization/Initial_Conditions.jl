module Initial_Conditions

using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Dyn_Data_Module
using ..Experiment_Configuration
using ..Variable_Mappings_Module
using ..Vert_Coordinate_Module
using ..Spectral_Dynamics_Module: Spectral_Initialize_Fields!, Get_Topography!
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

    elseif ic_name == :Moist_Spinup || ic_name == :Polvani_Test
        Init_3D_Standard!(mesh, atmo_data, dyn_data, vert_coord, config)

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
# 3. 3D Primitive Equation Initialization
#    (Logic extracted from exp/HSMoistt42/HS_Moist_Physics.jl)
# ==============================================================================
function Init_3D_Standard!(
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    vert_coord::Vert_Coordinate,
    config::Model_Config,
)
    # Default Reference Values
    sea_level_ps_ref = 1.0e5
    init_t = 264.0

    # 1. Set Topography (usually flat for HS, but function handles it)
    Get_Topography!(dyn_data.grid_geopots)

    # 2. Call Standard Spectral Initialization
    # This sets up T, lnps, and adds small perturbations to break symmetry
    Spectral_Initialize_Fields!(
        config,
        mesh,
        vert_coord,
        atmo_data,
        dyn_data,
        sea_level_ps_ref,
        init_t,
        dyn_data.grid_geopots,
    )
end

end
