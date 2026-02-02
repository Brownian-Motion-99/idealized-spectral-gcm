module Output_Manager_Module

using Base
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
    vert_coord=nothing
)
    if isfile(filename); rm(filename); end
    ds = NCDataset(filename, "c")

    # Dimensions
    defDim(ds, "lon",  mesh.nλ)
    defDim(ds, "lat",  mesh.nθ)
    defDim(ds, "time", Inf)

    if file_type == :plev
        defDim(ds, "plev", length(target_levels))
        v_lev = defVar(ds, "plev", Float64, ("plev",))
        v_lev[:] = target_levels
        v_lev.attrib["units"] = "Pa"
        
    elseif file_type == :raw
        # Only write Hybrid Coeffs if vert_coord is provided (3D case)
        defDim(ds, "lev", mesh.nd)
        
        if vert_coord !== nothing
            defDim(ds, "ilev", mesh.nd+1)
            
            # Calculate mid-point coefficients for metadata
            ak_m = 0.5 .* (vert_coord.ak[1:end-1] .+ vert_coord.ak[2:end])
            bk_m = 0.5 .* (vert_coord.bk[1:end-1] .+ vert_coord.bk[2:end])
            
            v_hyam = defVar(ds, "hyam", Float64, ("lev",))
            v_hyam[:] = ak_m
            v_hybm = defVar(ds, "hybm", Float64, ("lev",))
            v_hybm[:] = bk_m
        end
    end

    # Coords
    v_lon = defVar(ds, "lon", Float64, ("lon",))
    v_lon[:] = mesh.λc * 180/pi
    v_lat = defVar(ds, "lat", Float64, ("lat",))
    v_lat[:] = mesh.θc * 180/pi
    v_time = defVar(ds, "time", Float64, ("time",))

    # Variables
    for sym in requested_vars
        if !haskey(var_info_map, sym); continue; end
        nc_name, units, _, dim_code = var_info_map[sym]
        
        dims = if dim_code == 3
            (file_type == :plev) ? ("lon", "lat", "plev", "time") : ("lon", "lat", "lev", "time")
        else
            ("lon", "lat", "time")
        end
        
        defVar(ds, nc_name, Float64, dims, attrib=Dict("units"=>units))
    end

    return ds
end

# ==============================================================================
# 4. The Constructor
# ==============================================================================
function Output_Manager(
    mesh, 
    vert_coord, # Pass 'nothing' here for 2D runs in your driver script
    atmo_data,
    start_time, end_time,
    requested_vars;
    # Keywords
    filename::String="output.nc",
    do_raw_output::Bool=false,
    pressure_levels::Vector{Float64}=[100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    day_to_sec::Int64=86400,
    output_interval::Int64=86400,
    spinup_day::Float64=0.0,
    model_mode::Symbol=:PrimitiveEquation
)
    # A. Instantiate the Mode Type
    mode_obj = Mode_Factory(model_mode)

    nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
    spinup_time = start_time + Int64(spinup_day * day_to_sec)
    
    # Define Buffers (Initialize as empty by default)
    p3d_buf    = zeros(Float64, 0, 0, 0)
    interp_buf = zeros(Float64, 0, 0, 0)
    ak_m       = Float64[]
    bk_m       = Float64[]

    # B. 3D Specific Initialization
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

    # C. Initialize Files
    main_file_type = isa(mode_obj, PrimitiveEquationMode) ? :plev : :simple_2d
    
    # We pass vert_coord as a keyword now. If it's nothing (2D), it's ignored safely.
    ds_plev = Base.invokelatest(_Init_Single_File, filename, mesh, 
                                var_info_map, requested_vars, main_file_type, pressure_levels; 
                                vert_coord=vert_coord)

    ds_raw = nothing
    if do_raw_output && isa(mode_obj, PrimitiveEquationMode)
        raw_fn = replace(filename, ".nc" => "_raw.nc")
        ds_raw = Base.invokelatest(_Init_Single_File, raw_fn, mesh, 
                                   var_info_map, requested_vars, :raw; 
                                   vert_coord=vert_coord)
    end

    # D. Initialize Accumulators
    acc_main = Dict{Symbol, Array{Float64}}()
    
    for sym in requested_vars
        if !haskey(var_info_map, sym); continue; end
        _, _, _, dim_code = var_info_map[sym]
        
        if dim_code == 3
             # For 2D modes, this might create a 3D array if the user mistakenly requests a 3D var.
             # However, typically 2D modes only request 2D vars (dim_code=2).
             # If a 3D var is requested in 3D mode, we allocate (Nq, Nel, Nz).
             # If a 3D var is requested in 2D mode (e.g. theoretical), we still allocate it 
             # but the physics likely won't fill it.
            acc_main[sym] = zeros(Float64, nλ, nθ, nd)
        else
            acc_main[sym] = zeros(Float64, nλ, nθ)
        end
    end

    acc_raw  = deepcopy(acc_main)
    acc_plev = Dict{Symbol, Array{Float64}}() 

    return Output_Manager(
        nλ, nθ, nd, 
        atmo_data, mode_obj, 
        do_raw_output, pressure_levels, log.(pressure_levels),
        p3d_buf, interp_buf, ak_m, bk_m,
        day_to_sec, start_time, start_time, spinup_time, output_interval, 0, 1,
        ds_plev, ds_raw,
        acc_plev, acc_main, acc_raw, 
        requested_vars
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
    N     = Float64(mgr.sample_counter)
    var_info_map = Base.invokelatest(Get_Var_Info, Val(:PrimitiveEquation))

    # 1. Pre-calculate Averaged Pressure Grid (for Interpolation)
    ps_avg_2d = nothing
    if haskey(mgr.acc_main, :ps)
        ps_avg_2d = mgr.acc_main[:ps] ./ N
        Compute_Pressure_Grid!(mgr.p3d_buffer, mgr.vert_ak_mid, mgr.vert_bk_mid, ps_avg_2d)
    end

    # =========================================================================
    # 2. Process Main Output (Interpolated)
    # =========================================================================
    for sym in keys(mgr.acc_main)
        if !haskey(var_info_map, sym); continue; end
        nc_name = var_info_map[sym][1]
        
        # Robustly determine file variable dimensions
        nc_var = mgr.ds_plev[nc_name]
        file_ndim = ndims(nc_var) # e.g., 4 for (lon,lat,plev,time), 3 for (lon,lat,time)

        native_mean = mgr.acc_main[sym] ./ N
        
        # A. If File Variable is 4D (Pressure Levels + Time)
        if file_ndim == 4
            if ndims(native_mean) == 3
                t_ref = haskey(mgr.acc_main, :t) ? (mgr.acc_main[:t] ./ N) : nothing
                
                Interpolate_Field!(
                    mgr.interp_buffer, native_mean, mgr.p3d_buffer, ps_avg_2d,
                    mgr.log_targets, sym, mgr.atmo_data, t_ref
                )
                nc_var[:, :, :, t_idx] = mgr.interp_buffer
            else
                # Fallback: Trying to write 2D data to 4D var (rare warning case)
                @warn "Output_Manager: Dimension mismatch for $sym (File: 4D, Data: 2D). Skipping."
            end

        # B. If File Variable is 3D (Surface/2D + Time)
        elseif file_ndim == 3
            if ndims(native_mean) == 3
                # Flatten singleton if necessary (e.g. view(dat,:,:,1))
                nc_var[:, :, t_idx] = view(native_mean, :, :, 1)
            else
                nc_var[:, :, t_idx] = native_mean
            end
        end

        # Reset Accumulator
        fill!(mgr.acc_main[sym], 0.0)
    end
    NCDatasets.sync(mgr.ds_plev)

    # =========================================================================
    # 3. Process Raw Output (Native Grid)
    # =========================================================================
    if mgr.do_raw_output && mgr.ds_raw !== nothing
        for sym in keys(mgr.acc_raw)
            if !haskey(var_info_map, sym); continue; end # Safety check
            
            nc_name = var_info_map[sym][1]
            nc_var  = mgr.ds_raw[nc_name]
            
            # [FIX IS HERE] Check dimensionality before writing
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
        nc_name = var_info_map[sym][1]
        
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