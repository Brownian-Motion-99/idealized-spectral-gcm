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

# export Output_Manager, Update_Output!, Finalize_Output!
# export Lat_Lon_Pcolormesh, Zonal_Mean, Sigma_Zonal_Mean_Pcolormesh, Sigma_Zonal_Mean_Contourf

# mutable struct Output_Manager
#     nλ::Int64
#     nθ::Int64
#     nd::Int64
#     n_day::Int64
    
#     day_to_sec::Int64
#     start_time::Int64
#     end_time::Int64
#     current_time::Int64
#     spinup_day::Int64

#     λc::Array{Float64, 1}
#     θc::Array{Float64, 1}
#     σc::Array{Float64, 1}

#     n_daily_mean::Array{Float64, 1}    
#     ##########################################################################
#     # specral vor 
#     spe_vor_c_xyzt::Array{ComplexF64,4}
#     spe_vor_p_xyzt::Array{ComplexF64,4}
    
#     # specral div
#     spe_div_c_xyzt::Array{ComplexF64,4}
#     spe_div_p_xyzt::Array{ComplexF64,4}

#     # specral height or surface pressure
#     spe_lnps_c_xyzt::Array{ComplexF64,4}
#     spe_lnps_p_xyzt::Array{ComplexF64,4}

#     # specral temperature
#     spe_t_c_xyzt::Array{ComplexF64,4}
#     spe_t_p_xyzt::Array{ComplexF64,4}

#     # specral tracer
#     spe_tracers_c_xyzt::Array{ComplexF64,4}
#     spe_tracers_p_xyzt::Array{ComplexF64,4}
#     ##########################################################################
#     # grid w-e velocity
#     grid_u_n_xyzt::Array{Float64, 4}
#     grid_u_c_xyzt::Array{Float64, 4}
#     grid_u_p_xyzt::Array{Float64, 4}
    
#     # grid n-s velocity
#     grid_v_n_xyzt::Array{Float64, 4}
#     grid_v_c_xyzt::Array{Float64, 4}
#     grid_v_p_xyzt::Array{Float64, 4}

#     # grid surface pressure
#     grid_ps_c_xyzt::Array{Float64,4}
#     grid_ps_p_xyzt::Array{Float64,4}

#     # grid temperature
#     grid_t_n_xyzt::Array{Float64, 4}
#     grid_t_c_xyzt::Array{Float64, 4}
#     grid_t_p_xyzt::Array{Float64, 4}

#     # grid tracer
#     grid_tracers_n_xyzt::Array{Float64,4}
#     grid_tracers_c_xyzt::Array{Float64,4}
#     grid_tracers_p_xyzt::Array{Float64,4}

#     grid_tracers_diff_xyzt::Array{Float64,4}
#     grid_δtracers_xyzt::Array{Float64,4}
    
#     ##########################################################################
#     factor1_xyzt::Array{Float64,4}
#     factor2_xyzt::Array{Float64,4}
#     factor3_xyzt::Array{Float64,4}
#     factor4_xyzt::Array{Float64,4}

#     # pressure 
#     grid_p_full_xyzt::Array{Float64,4}
#     grid_p_half_xyzt::Array{Float64,4}
    
#     # geopotential
#     grid_geopots_xyzt::Array{Float64,4}
    
#     # Memory contrainer for temporal variables
#     # vor
#     grid_vor_xyzt::Array{Float64,4}
    
#     # div
#     grid_div_xyzt::Array{Float64,4}

#     # w-e velocity tendency
#     grid_δu_xyzt::Array{Float64,4}
    
#     # n-s velocity tendency
#     grid_δv_xyzt::Array{Float64,4}
#     #######################################################################
#     # equilibrium temperature in HS_Forcing
#     convection_xyzt::Array{Float64,4}

#     grid_z_full_xyzt::Array{Float64,4}
#     grid_w_full_xyzt::Array{Float64,4}
#     #######################################################################
#     grid_tracers_c_max_Tiffany_xyzt::Array{Float64,4}
#     grid_tracers_c_max_xyzt::Array{Float64,4}
    
    
    

    

# end

# function Output_Manager(mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, start_time::Int64, end_time::Int64, spinup_day::Int64, num_grid_tracters::Int64=1, num_spe_tracters::Int64=1) ### By CJY2
#     nλ = mesh.nλ
#     nθ = mesh.nθ
#     nd = mesh.nd
    
#     day_to_sec = 86400
#     current_time = start_time

#     λc = mesh.λc
#     θc = mesh.θc

#     #todo definition of sigma coordinate
#     bk = vert_coord.bk
#     σc = (bk[2:nd+1] + bk[1:nd])/2.0
  
#     n_day = Int64((end_time - start_time)/ (day_to_sec/4) )
#     n_daily_mean = zeros(Float64, n_day)

#     grid_geopots_xyzt = zeros(Float64, nλ, nθ, 1, n_day)
#     num_fourier, nθ, nd = 42, 64, 20
#     num_spherical = num_fourier + 1
#     #########################################################
#     # specral vor 
#     spe_vor_c_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     spe_vor_p_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)

#     # specral div
#     spe_div_c_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     spe_div_p_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)

#     # specral height or surface pressure
#     spe_lnps_c_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day)
#     spe_lnps_p_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, 1, n_day)

#     # specral temperature
#     spe_t_c_xyzt    = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     spe_t_p_xyzt    = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)

#     # specral tracer
#     spe_tracers_c_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     spe_tracers_p_xyzt = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     ##########################################################
#     # grid w-e velocity
#     grid_u_n_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
#     grid_u_c_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
#     grid_u_p_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)

#     # grid n-s velocity
#     grid_v_n_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
#     grid_v_c_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)
#     grid_v_p_xyzt  = zeros(Float64, nλ, nθ, nd, n_day)

#     # grid surface pressure
#     grid_ps_c_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)
#     grid_ps_p_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)

#     # grid temperature
#     grid_t_n_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_t_c_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_t_p_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)

#     # grid tracer
#     grid_tracers_n_xyzt     = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_tracers_c_xyzt     = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_tracers_p_xyzt     = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_tracers_diff_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_δtracers_xyzt      = zeros(Float64, nλ,  nθ, nd, n_day)
#     ############################################################
#     # factors
#     factor1_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 
#     factor2_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 
#     factor3_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 
#     factor4_xyzt = zeros(Float64, nλ,  nθ, nd, n_day) 

#     # pressure
#     grid_p_full_xyzt = zeros(Float64, nλ,  nθ, nd  , n_day) 
#     grid_p_half_xyzt = zeros(Float64, nλ,  nθ, nd+1, n_day) 
    
#     # geopotential
#     grid_geopots_xyzt = zeros(Float64, nλ,  nθ, 1, n_day)
#     ################################################
#     # spe_δvor_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     grid_vor_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day)
#     # grid_δvor_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    
#     # spe_δdiv_xyzt  = zeros(ComplexF64, num_fourier+1, num_spherical+1, nd, n_day)
#     grid_div_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day) 
#     # grid_δdiv_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    
#     # Tendency
#     grid_δu_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day) 
#     grid_δv_xyzt  = zeros(Float64, nλ,  nθ, nd, n_day) 

#     convection_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)

#     grid_z_full_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_w_full_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)

#     grid_tracers_c_max_Tiffany_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
#     grid_tracers_c_max_xyzt = zeros(Float64, nλ,  nθ, nd, n_day)
    
    
    
    
    
    
#     Output_Manager(nλ, nθ, nd, n_day,
#     day_to_sec, start_time, end_time, current_time, spinup_day,
#     λc, θc, σc, n_daily_mean, spe_vor_c_xyzt, spe_vor_p_xyzt, spe_div_c_xyzt, spe_div_p_xyzt, spe_lnps_c_xyzt, spe_lnps_p_xyzt, spe_t_c_xyzt, spe_t_p_xyzt, spe_tracers_c_xyzt, spe_tracers_p_xyzt, grid_u_n_xyzt, grid_u_c_xyzt, grid_u_p_xyzt, grid_v_n_xyzt, grid_v_c_xyzt, grid_v_p_xyzt, grid_ps_c_xyzt, grid_ps_p_xyzt, grid_t_n_xyzt, grid_t_c_xyzt, grid_t_p_xyzt, grid_tracers_n_xyzt, grid_tracers_c_xyzt, grid_tracers_p_xyzt, grid_tracers_diff_xyzt, grid_δtracers_xyzt, factor1_xyzt, factor2_xyzt, factor3_xyzt, factor4_xyzt, grid_p_full_xyzt, grid_p_half_xyzt, grid_geopots_xyzt, grid_vor_xyzt, grid_div_xyzt, grid_δu_xyzt, grid_δv_xyzt, convection_xyzt, grid_z_full_xyzt, grid_w_full_xyzt, grid_tracers_c_max_Tiffany_xyzt, grid_tracers_c_max_xyzt)
# end

# function Update_Output!(output_manager::Output_Manager, dyn_data::Dyn_Data, current_time::Int64)
#     @assert(current_time > output_manager.current_time)
#     output_manager.current_time = current_time
#     day_to_sec, start_time, n_day = output_manager.day_to_sec, output_manager.start_time, output_manager.n_day

#     n_daily_mean = output_manager.n_daily_mean
#     ############################################################
#     # specral vor 
#     spe_vor_c_xyzt  = output_manager.spe_vor_c_xyzt
#     spe_vor_p_xyzt  = output_manager.spe_vor_p_xyzt
    
#     # specral div
#     spe_div_c_xyzt  = output_manager.spe_div_c_xyzt
#     spe_div_p_xyzt  = output_manager.spe_div_p_xyzt

#     # specral height or surface pressure
#     spe_lnps_c_xyzt = output_manager.spe_lnps_c_xyzt
#     spe_lnps_p_xyzt = output_manager.spe_lnps_p_xyzt
    
#     # specral temperature
#     spe_t_c_xyzt = output_manager.spe_t_c_xyzt
#     spe_t_p_xyzt = output_manager.spe_t_p_xyzt
    
#     # specral tracer
#     spe_tracers_c_xyzt = output_manager.spe_tracers_c_xyzt
#     spe_tracers_p_xyzt = output_manager.spe_tracers_p_xyzt
#     ############################################################
#     # grid w-e velocity
#     grid_u_n_xyzt  = output_manager.grid_u_n_xyzt 
#     grid_u_c_xyzt  = output_manager.grid_u_c_xyzt 
#     grid_u_p_xyzt  = output_manager.grid_u_p_xyzt 
    
#     # grid n-s velocity
#     grid_v_n_xyzt  = output_manager.grid_v_n_xyzt
#     grid_v_c_xyzt  = output_manager.grid_v_c_xyzt
#     grid_v_p_xyzt  = output_manager.grid_v_p_xyzt

#     # grid surface pressure
#     grid_ps_c_xyzt = output_manager.grid_ps_c_xyzt
#     grid_ps_p_xyzt = output_manager.grid_ps_p_xyzt

#     # grid temperature
#     grid_t_n_xyzt  = output_manager.grid_t_n_xyzt
#     grid_t_c_xyzt  = output_manager.grid_t_c_xyzt
#     grid_t_p_xyzt  = output_manager.grid_t_p_xyzt

#     # grid tracer
#     grid_tracers_n_xyzt = output_manager.grid_tracers_n_xyzt
#     grid_tracers_c_xyzt = output_manager.grid_tracers_c_xyzt
#     grid_tracers_p_xyzt = output_manager.grid_tracers_p_xyzt

#     grid_tracers_diff_xyzt = output_manager.grid_tracers_diff_xyzt
#     grid_δtracers_xyzt     = output_manager.grid_δtracers_xyzt
#     ############################################################
#     # factors
#     factor1_xyzt = output_manager.factor1_xyzt
#     factor2_xyzt = output_manager.factor2_xyzt
#     factor3_xyzt = output_manager.factor3_xyzt
#     factor4_xyzt = output_manager.factor4_xyzt

    
#     # pressure
#     grid_p_full_xyzt = output_manager.grid_p_full_xyzt
#     grid_p_half_xyzt = output_manager.grid_p_half_xyzt
    
#     # geopotential
#     grid_geopots_xyzt = output_manager.grid_geopots_xyzt
#     ############################################################
#     grid_vor_xyzt  = output_manager.grid_vor_xyzt
    
#     grid_div_xyzt  = output_manager.grid_div_xyzt
    
#     # Tendency
#     grid_δu_xyzt   = output_manager.grid_δu_xyzt
#     grid_δv_xyzt   = output_manager.grid_δv_xyzt

#     convection_xyzt = output_manager.convection_xyzt

#     grid_z_full_xyzt = output_manager.grid_z_full_xyzt
#     grid_w_full_xyzt = output_manager.grid_w_full_xyzt

#     grid_tracers_c_max_Tiffany_xyzt = output_manager.grid_tracers_c_max_Tiffany_xyzt
#     grid_tracers_c_max_xyzt = output_manager.grid_tracers_c_max_xyzt
    
#     i_day = Int(div(current_time - start_time - 1, day_to_sec/4) + 1)

#     if(i_day > n_day)
#         @info "Warning: i_day > n_day in Output_Manager:Update!"
#         return 
#     end
    
#     ############################################################
#     # specral vor 
#     spe_vor_c_xyzt[:,:,:,i_day] .= dyn_data.spe_vor_c[:,:,:]
#     spe_vor_p_xyzt[:,:,:,i_day] .= dyn_data.spe_vor_p[:,:,:]
        
#     # specral div
#     spe_div_c_xyzt[:,:,:,i_day] .= dyn_data.spe_div_c[:,:,:]
#     spe_div_p_xyzt[:,:,:,i_day] .= dyn_data.spe_div_p[:,:,:]
    
#     # specral height or surface pressure
#     spe_lnps_c_xyzt[:,:,:,i_day] .= dyn_data.spe_lnps_c[:,:,:]
#     spe_lnps_p_xyzt[:,:,:,i_day] .= dyn_data.spe_lnps_p[:,:,:]
    
#     # specral temperature
#     spe_t_c_xyzt[:,:,:,i_day] .= dyn_data.spe_t_c[:,:,:]
#     spe_t_p_xyzt[:,:,:,i_day] .= dyn_data.spe_t_p[:,:,:]
        
#     # specral tracer
#     spe_tracers_c_xyzt[:,:,:,i_day] .= dyn_data.spe_tracers_c[:,:,:]
#     spe_tracers_p_xyzt[:,:,:,i_day] .= dyn_data.spe_tracers_p[:,:,:]
#     ############################################################
#     # grid w-e velocity
#     grid_u_n_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_n[:,:,:] 
#     grid_u_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_c[:,:,:] 
#     grid_u_p_xyzt[:,:,:,i_day]  .= dyn_data.grid_u_p[:,:,:]
    
#     # grid n-s velocity
#     grid_v_n_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_n[:,:,:]
#     grid_v_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_c[:,:,:]
#     grid_v_p_xyzt[:,:,:,i_day]  .= dyn_data.grid_v_p[:,:,:]

#     # grid surface pressure
#     grid_ps_c_xyzt[:,:,:,i_day] .= dyn_data.grid_ps_c[:,:,:]
#     grid_ps_p_xyzt[:,:,:,i_day] .= dyn_data.grid_ps_p[:,:,:]

#     # grid temperature
#     grid_t_n_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_n[:,:,:]
#     grid_t_c_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_c[:,:,:]
#     grid_t_p_xyzt[:,:,:,i_day]  .= dyn_data.grid_t_p[:,:,:]

#     # grid tracer
#     grid_tracers_n_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_n[:,:,:]
#     grid_tracers_c_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_c[:,:,:]
#     grid_tracers_p_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_p[:,:,:]

#     grid_tracers_diff_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_diff[:,:,:]
#     grid_δtracers_xyzt[:,:,:,i_day]     .= dyn_data.grid_δtracers[:,:,:]
#     ############################################################
#     # factors
#     factor1_xyzt[:,:,:,i_day] .= dyn_data.factor1[:,:,:]
#     factor2_xyzt[:,:,:,i_day] .= dyn_data.factor2[:,:,:]
#     factor3_xyzt[:,:,:,i_day] .= dyn_data.factor3[:,:,:]
#     factor4_xyzt[:,:,:,i_day] .= dyn_data.factor4[:,:,:]

    
#     # pressure
#     grid_p_full_xyzt[:,:,:,i_day] .= dyn_data.grid_p_full[:,:,:]
#     grid_p_half_xyzt[:,:,:,i_day] .= dyn_data.grid_p_half[:,:,:]
    
#     # geopotential
#     grid_geopots_xyzt[:,:,:,i_day] .= dyn_data.grid_geopots[:,:,:]
#     ############################################################
#     grid_vor_xyzt[:,:,:,i_day]  .= dyn_data.grid_vor[:,:,:]
#     grid_div_xyzt[:,:,:,i_day]  .= dyn_data.grid_div[:,:,:]
    
#     # Tendency
#     grid_δu_xyzt[:,:,:,i_day]   .= dyn_data.grid_δu[:,:,:]
#     grid_δv_xyzt[:,:,:,i_day]   .= dyn_data.grid_δv[:,:,:]

#     convection_xyzt[:,:,:,i_day] .= dyn_data.convection[:,:,:]

#     grid_z_full_xyzt[:,:,:,i_day] .= dyn_data.grid_z_full[:,:,:]
#     grid_w_full_xyzt[:,:,:,i_day] .= dyn_data.grid_w_full[:,:,:]

#     grid_tracers_c_max_Tiffany_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_c_max_Tiffany[:,:,:]
#     grid_tracers_c_max_xyzt[:,:,:,i_day] .= dyn_data.grid_tracers_c_max[:,:,:]
    

#     n_daily_mean[i_day] += 1
# end

# function Finalize_Output!(output_manager::Output_Manager, save_file_name::String = "None", mean_save_file_name::String = "None")

#     n_day = output_manager.n_day

#     ############################################################    
#     # specral vor 
#     spe_vor_c_xyzt  = output_manager.spe_vor_c_xyzt
#     spe_vor_p_xyzt  = output_manager.spe_vor_p_xyzt
    
#     # specral div
#     spe_div_c_xyzt  = output_manager.spe_div_c_xyzt
#     spe_div_p_xyzt  = output_manager.spe_div_p_xyzt

#     # specral height or surface pressure
#     spe_lnps_c_xyzt = output_manager.spe_lnps_c_xyzt
#     spe_lnps_p_xyzt = output_manager.spe_lnps_p_xyzt
    
#     # specral temperature
#     spe_t_c_xyzt    = output_manager.spe_t_c_xyzt
#     spe_t_p_xyzt    = output_manager.spe_t_p_xyzt
    
#     # specral tracer
#     spe_tracers_c_xyzt = output_manager.spe_tracers_c_xyzt
#     spe_tracers_p_xyzt = output_manager.spe_tracers_p_xyzt
#     ############################################################
#     # grid w-e velocity
#     grid_u_n_xyzt  = output_manager.grid_u_n_xyzt 
#     grid_u_c_xyzt  = output_manager.grid_u_c_xyzt 
#     grid_u_p_xyzt  = output_manager.grid_u_p_xyzt 
    
#     # grid n-s velocity
#     grid_v_n_xyzt  = output_manager.grid_v_n_xyzt
#     grid_v_c_xyzt  = output_manager.grid_v_c_xyzt
#     grid_v_p_xyzt  = output_manager.grid_v_p_xyzt

#     # grid surface pressure
#     grid_ps_c_xyzt = output_manager.grid_ps_c_xyzt
#     grid_ps_p_xyzt = output_manager.grid_ps_p_xyzt

#     # grid temperature
#     grid_t_n_xyzt  = output_manager.grid_t_n_xyzt
#     grid_t_c_xyzt  = output_manager.grid_t_c_xyzt
#     grid_t_p_xyzt  = output_manager.grid_t_p_xyzt

#     # grid tracer
#     grid_tracers_n_xyzt = output_manager.grid_tracers_n_xyzt  
#     grid_tracers_c_xyzt = output_manager.grid_tracers_c_xyzt
#     grid_tracers_p_xyzt = output_manager.grid_tracers_p_xyzt

#     grid_tracers_diff_xyzt = output_manager.grid_tracers_diff_xyzt
#     grid_δtracers_xyzt     = output_manager.grid_δtracers_xyzt
#     ############################################################
#     # factors
#     factor1_xyzt = output_manager.factor1_xyzt
#     factor2_xyzt = output_manager.factor2_xyzt
#     factor3_xyzt = output_manager.factor3_xyzt
#     factor4_xyzt = output_manager.factor4_xyzt


#     # pressure
#     grid_p_full_xyzt = output_manager.grid_p_full_xyzt
#     grid_p_half_xyzt = output_manager.grid_p_half_xyzt
    
#     # geopotential
#     grid_geopots_xyzt = output_manager.grid_geopots_xyzt
#     ############################################################
#     grid_vor_xyzt  = output_manager.grid_vor_xyzt    
#     grid_div_xyzt  = output_manager.grid_div_xyzt
    
#     # Tendency
#     grid_δu_xyzt   = output_manager.grid_δu_xyzt
#     grid_δv_xyzt   = output_manager.grid_δv_xyzt

#     convection_xyzt = output_manager.convection_xyzt

#     grid_z_full_xyzt = output_manager.grid_z_full_xyzt
#     grid_w_full_xyzt = output_manager.grid_w_full_xyzt
#     grid_tracers_c_max_Tiffany_xyzt = output_manager.grid_tracers_c_max_Tiffany_xyzt
#     grid_tracers_c_max_xyzt = output_manager.grid_tracers_c_max_xyzt
    
#     ##############################################################################
#     spe_vor_c_final = spe_vor_c_xyzt[:,:,:,end]
#     spe_vor_p_final = spe_vor_p_xyzt[:,:,:,end]
#     spe_div_c_final = spe_div_c_xyzt[:,:,:,end]
#     spe_div_p_final = spe_div_p_xyzt[:,:,:,end]

#     spe_lnps_c_final = spe_lnps_c_xyzt[:,:,:,end]
#     spe_lnps_p_final = spe_lnps_p_xyzt[:,:,:,end]

#     spe_t_c_final = spe_t_c_xyzt[:,:,:,end]
#     spe_t_p_final = spe_t_p_xyzt[:,:,:,end]

#     spe_tracers_c_final = spe_tracers_c_xyzt[:,:,:,end]
#     spe_tracers_p_final = spe_tracers_p_xyzt[:,:,:,end]

#     grid_tracers_n_final = grid_tracers_n_xyzt[:,:,:,end]

#     grid_u_n_final = grid_u_n_xyzt[:,:,:,end]
#     grid_u_c_final = grid_u_c_xyzt[:,:,:,end]
#     grid_u_p_final = grid_u_p_xyzt[:,:,:,end]

#     grid_v_n_final = grid_v_n_xyzt[:,:,:,end]
#     grid_v_c_final = grid_v_c_xyzt[:,:,:,end]
#     grid_v_p_final = grid_v_p_xyzt[:,:,:,end]

#     grid_ps_c_final = grid_ps_c_xyzt[:,:,:,end]
#     grid_ps_p_final = grid_ps_p_xyzt[:,:,:,end]

#     grid_t_n_final = grid_t_n_xyzt[:,:,:,end]
#     grid_t_c_final = grid_t_c_xyzt[:,:,:,end]
#     grid_t_p_final = grid_t_p_xyzt[:,:,:,end]

#     grid_tracers_c_final = grid_tracers_c_xyzt[:,:,:,end]
#     grid_tracers_p_final = grid_tracers_p_xyzt[:,:,:,end]

#     grid_tracers_diff_final = grid_tracers_diff_xyzt[:,:,:,end]
#     grid_δtracers_final = grid_δtracers_xyzt[:,:,:,end]

#     grid_p_full_final = grid_p_full_xyzt[:,:,:,end]
#     grid_p_half_final = grid_p_half_xyzt[:,:,:,end]

#     grid_geopots_final = grid_geopots_xyzt[:,:,:,end]

#     grid_vor_final = grid_vor_xyzt[:,:,:,end]
#     grid_div_final = grid_div_xyzt[:,:,:,end]

#     grid_δu_final = grid_δu_xyzt[:,:,:,end]
#     grid_δv_final = grid_δv_xyzt[:,:,:,end]

#     convection_final = convection_xyzt[:,:,:,end]

#     grid_w_full_final = grid_w_full_xyzt[:,:,:,end]
    
    
    
    
    
    
    
#     if save_file_name != "None"
#         @save save_file_name spe_vor_c_final spe_vor_p_final spe_div_c_final spe_div_p_final spe_lnps_c_final spe_lnps_p_final spe_t_c_final spe_t_p_final spe_tracers_c_final spe_tracers_p_final grid_tracers_n_final grid_u_n_final grid_u_c_final grid_u_p_final grid_v_n_final grid_v_c_final grid_v_p_final grid_ps_c_final grid_ps_p_final grid_t_n_final grid_t_c_final grid_t_p_final grid_tracers_c_final grid_tracers_p_final grid_tracers_diff_final grid_δtracers_final grid_p_full_final grid_p_half_final grid_geopots_final grid_vor_final grid_div_final grid_δu_final grid_δv_final convection_final grid_w_full_final 
#     end

#     if mean_save_file_name != "None"
#         @save mean_save_file_name spe_vor_c_xyzt spe_vor_p_xyzt spe_div_c_xyzt spe_div_p_xyzt spe_lnps_c_xyzt spe_lnps_p_xyzt spe_t_c_xyzt spe_t_p_xyzt spe_tracers_c_xyzt spe_tracers_p_xyzt grid_tracers_n_xyzt grid_u_n_xyzt grid_u_c_xyzt  grid_u_p_xyzt grid_v_n_xyzt grid_v_c_xyzt grid_v_p_xyzt grid_ps_c_xyzt grid_ps_p_xyzt grid_t_n_xyzt grid_t_c_xyzt grid_t_p_xyzt grid_tracers_c_xyzt grid_tracers_p_xyzt grid_tracers_diff_xyzt grid_δtracers_xyzt factor1_xyzt factor2_xyzt factor3_xyzt grid_p_full_xyzt grid_p_half_xyzt grid_geopots_xyzt  grid_vor_xyzt  grid_div_xyzt grid_δu_xyzt grid_δv_xyzt convection_xyzt factor4_xyzt grid_z_full_xyzt grid_w_full_xyzt grid_tracers_c_max_Tiffany_xyzt grid_tracers_c_max_xyzt
#     end
# end