module Output_Manager_Module

using Infiltrator
using Base
using Dates
using NCDatasets
using Statistics

using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Dyn_Data_Module
using ..Vert_Coordinate_Module
using ..Vertical_Interpolation_Module
using ..Output_Mappings_Module
using ..Variable_Mappings_Module

export Output_Manager, Update_Output!, Finalize_Output!
export PrimitiveEquationMode, BarotropicMode, ShallowWaterMode

# ==============================================================================
# 1. Type Hierarchy (The Backbone of Dispatch)
# ==============================================================================
abstract type AbstractModelMode end

struct PrimitiveEquationMode <: AbstractModelMode end
struct BarotropicMode <: AbstractModelMode end
struct ShallowWaterMode <: AbstractModelMode end

# Helper to map legacy symbols to Types
function Mode_Factory(sym::Symbol)
    if sym == :PrimitiveEquation
        return PrimitiveEquationMode()
    elseif sym == :Barotropic
        return BarotropicMode()
    elseif sym == :ShallowWater
        return ShallowWaterMode()
    else
        error("Output_Manager: Unknown model_mode symbol :$sym")
    end
end

# Helper to map Types back to Symbols (for file attributes/legacy APIs)
mode_symbol(::PrimitiveEquationMode) = :PrimitiveEquation
mode_symbol(::BarotropicMode)        = :Barotropic
mode_symbol(::ShallowWaterMode)      = :ShallowWater

# ==============================================================================
# 2. The Parameterized Struct
# ==============================================================================
mutable struct Output_Manager{M <: AbstractModelMode}
    # Grid Metadata
    nλ::Int64
    nθ::Int64
    nd::Int64

    # --- Physics Context ---
    atmo_data::Atmo_Data
    mode::M  # <--- The concrete type (PrimitiveEquationMode, etc.) replaces the Symbol

    # --- Configuration ---
    do_raw_output::Bool              
    target_levels::Vector{Float64}   
    log_targets::Vector{Float64}     
    
    # --- Buffers (Pre-allocated) ---
    p3d_buffer::Array{Float64, 3}    
    interp_buffer::Array{Float64, 3} 
    
    vert_ak_mid::Vector{Float64}
    vert_bk_mid::Vector{Float64}

    # --- Time Control ---
    day_to_sec::Int64
    start_time::Int64
    current_time::Int64
    spinup_time::Int64
    output_interval::Int64
    sample_counter::Int64
    output_index::Int64

    # --- File Handles ---
    ds_plev::Union{NCDataset, Nothing} 
    ds_raw::Union{NCDataset, Nothing}  

    # --- Accumulators ---
    # Path A: We accumulate on the native grid in acc_main
    acc_plev::Dict{Symbol, Array{Float64}}       
    acc_main::Dict{Symbol, Array{Float64}}
    acc_raw::Dict{Symbol, Array{Float64}}        
    
    active_symbols::Vector{Symbol}
end

# ==============================================================================
# 3. Internal Initialization Logic
# ==============================================================================
function _Init_Single_File(
    filename, mesh, var_info_map, 
    requested_vars, file_type::Symbol, target_levels=Float64[];
    vert_coord=nothing,
    global_meta::Dict{String, String}=Dict{String, String}()
)
    if isfile(filename); rm(filename); end
    ds = NCDataset(filename, "c")

    # --- Global Attributes (CF-1.11 Compliance) --- #
    for (key, val) in global_meta
        ds.attrib[key] = val
    end
    ds.attrib["Conventions"] = "CF-1.11"
    ds.attrib["history"]     = "Created $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))"
    ds.attrib["source"]      = "idealized-spectral-gcm"
    # --- Global Attributes (CF-1.11 Compliance) --- #

    # --- Dimensions --- #
    defDim(ds, "lon",  mesh.nλ)
    defDim(ds, "lat",  mesh.nθ)
    defDim(ds, "time", Inf)
    # --- Dimensions --- #

    # --- Coordinate Variables & Vertical Metadata --- #
    # Longitude
    v_lon = defVar(ds, "lon", Float64, ("lon",), 
                   attrib=Dict("units"=>"degrees_east", "standard_name"=>"longitude"))
    v_lon[:] = mesh.λc * 180/pi

    # Latitude
    v_lat = defVar(ds, "lat", Float64, ("lat",), 
                   attrib=Dict("units"=>"degrees_north", "standard_name"=>"latitude"))
    v_lat[:] = mesh.θc * 180/pi

    # Time
    v_time = defVar(ds, "time", Float64, ("time",), 
                    attrib=Dict("units"=>"days since 2000-01-01 00:00:00", "calendar"=>"proleptic_gregorian"))

    if file_type == :plev
        defDim(ds, "plev", length(target_levels))
        v_lev = defVar(ds, "plev", Float64, ("plev",), 
                       attrib=Dict("units"=>"Pa", "standard_name"=>"air_pressure", "positive"=>"down"))
        v_lev[:] = target_levels
        
    elseif file_type == :raw && vert_coord !== nothing
        # Hybrid-Sigma coordinate definition
        defDim(ds, "lev", mesh.nd)
        v_lev = defVar(ds, "lev", Float64, ("lev",))
        v_lev[:] = collect(1:mesh.nd)
        
        # Attributes that tell post-processors how to compute 3D Pressure
        v_lev.attrib["standard_name"] = "atmosphere_hybrid_sigma_pressure_coordinate"
        v_lev.attrib["formula_terms"] = "ap: hyam b: hybm ps: ps p0: p0"
        v_lev.attrib["positive"]      = "down"

        # Reference Pressure
        v_p0 = defVar(ds, "p0", Float64, (), attrib=Dict("units"=>"Pa"))
        v_p0[:] = vert_coord.p_ref

        # Interface coefficients (Mid-points for data layers)
        ak_m = 0.5 .* (vert_coord.ak[1:end-1] .+ vert_coord.ak[2:end])
        bk_m = 0.5 .* (vert_coord.bk[1:end-1] .+ vert_coord.bk[2:end])

        defVar(ds, "hyam", Float64, ("lev",), attrib=Dict("units"=>"Pa", "long_name"=>"hybrid A coefficient at layer midpoints"))[1:end] = ak_m
        defVar(ds, "hybm", Float64, ("lev",), attrib=Dict("units"=>"1",  "long_name"=>"hybrid B coefficient at layer midpoints"))[1:end] = bk_m
    end
    # --- Coordinate Variables & Vertical Metadata --- #

    # --- Data Variables --- #
    for sym in requested_vars
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym] # Access via VarMeta struct
        
        dims = if meta.dims == 3
            (file_type == :plev) ? ("lon", "lat", "plev", "time") : ("lon", "lat", "lev", "time")
        else
            ("lon", "lat", "time")
        end
        
        v = defVar(ds, meta.nc_name, Float64, dims)
        v.attrib["units"]         = meta.units
        v.attrib["long_name"]     = meta.long_name
        v.attrib["standard_name"] = meta.std_name
    end
    # --- Data Variables --- #

    return ds
end

# ==============================================================================
# 4. The Constructor
# ==============================================================================
function Output_Manager(
    mesh, vert_coord, atmo_data, 
    start_time, end_time, requested_vars;
    filename::String = "output.nc",
    do_raw_output::Bool = true,
    pressure_levels::Vector{Float64} = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    day_to_sec::Int64 = 86400,
    output_interval::Int64 = 86400,
    spinup_day::Float64 = 0.0,
    model_mode::Symbol = :PrimitiveEquation,
    institute::String = "My Research Lab",
    experiment_id::String = "Simulation_v1"
)
    # --- Instantiate the Mode Type --- #
    mode_obj    = Mode_Factory(model_mode)
    nλ, nθ, nd  = mesh.nλ, mesh.nθ, mesh.nd
    spinup_time = start_time + Int64(spinup_day * day_to_sec)

    # Global Metadata Dictionary for the NetCDF header
    global_meta = Dict(
        "title"       => "Output of idealized-spectral-gcm",
        "institution" => institute,
        "experiment"  => experiment_id,
        "references"  => "https://github.com/Brownian-Motion-99/idealized-spectral-gcm.git"
    )
    
    # Define Buffers (Initialize as empty by default)
    p3d_buf    = zeros(Float64, 0, 0, 0)
    interp_buf = zeros(Float64, 0, 0, 0)
    ak_m       = Float64[]
    bk_m       = Float64[]
    # --- Instantiate the Mode Type --- #

    # --- 3D Specific Initialization --- #
    if isa(mode_obj, PrimitiveEquationMode)
        # Validation
        if !(:ps in requested_vars); push!(requested_vars, :ps); end
        if (:z in requested_vars) && !(:t in requested_vars); error("Need :t for :z"); end
        
        if vert_coord === nothing
             error("Output_Manager: vert_coord is required for :PrimitiveEquation mode")
        end

        # Allocation (Only happens for 3D)
        n_plev     = length(pressure_levels)
        p3d_buf    = zeros(Float64, nλ, nθ, nd)
        interp_buf = zeros(Float64, nλ, nθ, n_plev)
        
        ak_m = 0.5 .* (vert_coord.ak[1:end-1] .+ vert_coord.ak[2:end])
        bk_m = 0.5 .* (vert_coord.bk[1:end-1] .+ vert_coord.bk[2:end])
    end

    # Fetch variable info
    var_info_map = Base.invokelatest(Get_Var_Info, Val(mode_symbol(mode_obj)))
    # --- 3D Specific Initialization --- #

    # --- Initialize Files --- #
    main_file_type = isa(mode_obj, PrimitiveEquationMode) ? :plev : :simple_2d
    ds_plev = Base.invokelatest(
        _Init_Single_File, 
        filename, mesh, 
        var_info_map, requested_vars, main_file_type, pressure_levels; 
        vert_coord=vert_coord, global_meta=global_meta
    )

    ds_raw = nothing
    if do_raw_output && isa(mode_obj, PrimitiveEquationMode)
        raw_fn = replace(filename, ".nc" => "_raw.nc")
        ds_raw = Base.invokelatest(
            _Init_Single_File, 
            raw_fn, mesh, 
            var_info_map, requested_vars, :raw; 
            vert_coord=vert_coord, global_meta=global_meta
        )
    end
    # --- Initialize Files --- #

    # --- Initialize Accumulators --- #
    acc_main = Dict{Symbol, Array{Float64}}()
    for sym in requested_vars
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym]
        acc_main[sym] = (meta.dims == 3) ? zeros(Float64, nλ, nθ, nd) : zeros(Float64, nλ, nθ)
    end

    acc_raw  = deepcopy(acc_main)
    acc_plev = Dict{Symbol, Array{Float64}}() 
    # --- Initialize Accumulators --- #

    return Output_Manager(
        nλ, nθ, nd, atmo_data, mode_obj, 
        do_raw_output, pressure_levels, log.(pressure_levels),
        p3d_buf, interp_buf, ak_m, bk_m,
        day_to_sec, start_time, start_time, spinup_time, output_interval, 0, 1,
        ds_plev, ds_raw, acc_plev, acc_main, acc_raw, requested_vars
    )
end

# ==============================================================================
# 5. Runtime Logic (Dispatch Enabled)
# ==============================================================================

# Generic Wrapper
function Update_Output!(manager::Output_Manager{M}, dyn_data, current_time::Int64) where M
    manager.current_time = current_time
    if current_time <= manager.spinup_time; return; end

    # Fetch Data using the Mode Type (cleaner dispatch)
    # Assuming Get_Dyn_Var_Map still expects Val{Symbol} for now:
    live_data = Get_Dyn_Var_Map(dyn_data, Val(mode_symbol(manager.mode)))
    
    # CALL SPECIALIZED ACCUMULATION
    _accumulate_core!(manager.mode, manager, live_data)

    manager.sample_counter += 1

    time_elapsed = current_time - manager.start_time
    if time_elapsed > 0 && (time_elapsed % manager.output_interval == 0)
        Flush_to_Disk!(manager)
    end
end

# --- Specialized Accumulation: Primitive Equation ---
function _accumulate_core!(::PrimitiveEquationMode, mgr, live_data)
    for sym in mgr.active_symbols
        if !haskey(live_data, sym); continue; end
        
        raw_data = live_data[sym]
        
        # Main Accumulator
        if haskey(mgr.acc_main, sym)
            _accumulate_buffer!(mgr.acc_main[sym], raw_data)
        end

        # Raw Accumulator (if distinct file requested)
        if mgr.do_raw_output && haskey(mgr.acc_raw, sym)
            _accumulate_buffer!(mgr.acc_raw[sym], raw_data)
        end
    end
end

# --- Specialized Accumulation: 2D Modes ---
function _accumulate_core!(::AbstractModelMode, mgr, live_data)
    # Generic fallback for Barotropic/ShallowWater
    for sym in mgr.active_symbols
        if !haskey(live_data, sym); continue; end
        
        raw_data = live_data[sym]
        
        if haskey(mgr.acc_main, sym)
            # 2D models usually don't have dimension mismatches, direct add
            mgr.acc_main[sym] .+= raw_data
        end
    end
end

# Helper to handle singleton dimensions (e.g. Surface Pressure in 3D array)
@inline function _accumulate_buffer!(buffer, data)
    if ndims(buffer) == ndims(data)
        buffer .+= data
    elseif ndims(buffer) == 2 && ndims(data) == 3
        buffer .+= view(data, :, :, 1)
    end
end

# ==============================================================================
# 6. IO Logic (Dispatch Enabled)
# ==============================================================================

function Flush_to_Disk!(manager::Output_Manager{M}) where M
    if manager.sample_counter == 0; return; end

    manager.ds_plev["time"][manager.output_index] = manager.current_time / manager.day_to_sec
    
    if manager.do_raw_output && manager.ds_raw !== nothing
        manager.ds_raw["time"][manager.output_index] = manager.current_time / manager.day_to_sec
    end
    
    # Dispatch to specific writing logic
    _write_core!(manager.mode, manager)

    manager.sample_counter = 0
    manager.output_index += 1
end

# --- Writer: Primitive Equation (Includes Interpolation & Robust Dimension Check) ---
function _write_core!(::PrimitiveEquationMode, mgr)
    t_idx = mgr.output_index
    N            = Float64(mgr.sample_counter)
    var_info_map = Base.invokelatest(Get_Var_Info, Val(:PrimitiveEquation))

    # Calculate Interpolation Grid
    ps_avg_2d = mgr.acc_main[:ps] ./ N
    Compute_Pressure_Grid!(mgr.p3d_buffer, mgr.vert_ak_mid, mgr.vert_bk_mid, ps_avg_2d)

    for sym in keys(mgr.acc_main)
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym]
        nc_var = mgr.ds_plev[meta.nc_name]
        native_mean = mgr.acc_main[sym] ./ N
        
        if ndims(nc_var) == 4 # (lon, lat, plev, time)
            t_ref = haskey(mgr.acc_main, :t) ? (mgr.acc_main[:t] ./ N) : nothing
            Interpolate_Field!(mgr.interp_buffer, native_mean, mgr.p3d_buffer, ps_avg_2d,
                               mgr.log_targets, sym, mgr.atmo_data, t_ref)
            nc_var[:, :, :, t_idx] = mgr.interp_buffer
        else
            nc_var[:, :, t_idx] = (ndims(native_mean) == 3) ? view(native_mean, :, :, 1) : native_mean
        end
        fill!(mgr.acc_main[sym], 0.0)
    end
    NCDatasets.sync(mgr.ds_plev)

    # =========================================================================
    # 3. Process Raw Output (Native Grid)
    # =========================================================================
    if mgr.do_raw_output && mgr.ds_raw !== nothing
        for sym in keys(mgr.acc_raw)
            if !haskey(var_info_map, sym); continue; end # Safety check
            
            @infiltrate
            nc_name = var_info_map[sym].nc_name
            nc_var  = mgr.ds_raw[nc_name]
            
            # Check dimensionality before writing
            file_ndim = ndims(nc_var) 
            raw_mean  = mgr.acc_raw[sym] ./ N

            if file_ndim == 4 # (Lon, Lat, Lev, Time)
                nc_var[:, :, :, t_idx] = raw_mean
            elseif file_ndim == 3 # (Lon, Lat, Time) - e.g. Surface Pressure
                if ndims(raw_mean) == 3
                    nc_var[:, :, t_idx] = view(raw_mean, :, :, 1)
                else
                    nc_var[:, :, t_idx] = raw_mean
                end
            end
            
            fill!(mgr.acc_raw[sym], 0.0)
        end
        NCDatasets.sync(mgr.ds_raw)
    end
end

# --- Writer: 2D Modes (Direct Dump) ---
function _write_core!(::AbstractModelMode, mgr)
    t_idx = mgr.output_index
    N     = Float64(mgr.sample_counter)
    # We use mode_symbol to get the generic map for Barotropic/ShallowWater
    var_info_map = Base.invokelatest(Get_Var_Info, Val(mode_symbol(mgr.mode)))

    for sym in keys(mgr.acc_main)
        if !haskey(var_info_map, sym); continue; end
        nc_name = var_info_map[sym].nc_name
        
        # Direct write, no interpolation
        mgr.ds_plev[nc_name][:, :, t_idx] = mgr.acc_main[sym] ./ N
        fill!(mgr.acc_main[sym], 0.0)
    end
    NCDatasets.sync(mgr.ds_plev)
end

function Finalize_Output!(manager::Output_Manager)
    if manager.ds_plev !== nothing; close(manager.ds_plev); end
    if manager.ds_raw !== nothing; close(manager.ds_raw); end
end

end