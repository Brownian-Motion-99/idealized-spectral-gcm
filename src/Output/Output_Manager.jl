module Output_Manager_Module

using Base
using Dates
using NCDatasets

using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Dyn_Data_Module
using ..Vert_Coordinate_Module
using ..Vertical_Interpolation_Module
using ..Output_Mappings_Module
using ..Variable_Mappings_Module

export Output_Manager, Update_Output!, Finalize_Output!, Rotate_NC_Chunk!
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
    mode::M

    # --- Configuration ---
    do_plev_output::Bool             # true = also write pressure-level interpolated output
    target_levels::Vector{Float64}
    log_targets::Vector{Float64}

    # --- Buffers (Pre-allocated, 3D only) ---
    p3d_buffer::Array{Float64, 3}
    interp_buffer::Array{Float64, 3}

    vert_ak::Vector{Float64}
    vert_bk::Vector{Float64}

    # --- Time Control ---
    day_to_sec::Int64
    start_time::Int64
    current_time::Int64
    spinup_time::Int64
    output_interval::Int64
    sample_counter::Int64
    output_index::Int64

    # --- File Handles ---
    ds_plev::Union{NCDataset, Nothing}  # pressure-level output (optional, 3D only)
    ds_raw::Union{NCDataset, Nothing}   # native sigma/hybrid (or 2D native) — always created

    # --- Accumulators ---
    # acc_plev: pressure-level accumulator (only allocated when do_plev_output = true, 3D only)
    acc_plev::Dict{Symbol, Array{Float64}}
    # acc_raw:  native grid accumulator (always allocated, used for all model types)
    acc_raw::Dict{Symbol, Array{Float64}}

    active_symbols::Vector{Symbol}
    var_info_map::Dict{Symbol, VarMeta}

    # --- Chunk I/O State ---
    base_filename::String            # base path without extension, e.g. "exp/HSt42/output"
    mesh_λc::Vector{Float64}        # longitude grid in radians (for file re-creation on rotation)
    mesh_θc::Vector{Float64}        # latitude grid in radians
    vert_coord_cache::Any           # Vert_Coordinate or nothing
    global_meta::Dict{String, String}  # NC global attributes
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
    defVar(ds, "time", Float64, ("time",),
           attrib=Dict("units"=>"days since 2000-01-01 00:00:00", "calendar"=>"proleptic_gregorian"))

    if file_type == :plev
        defDim(ds, "plev", length(target_levels))
        v_lev = defVar(ds, "plev", Float64, ("plev",),
                       attrib=Dict("units"=>"Pa", "standard_name"=>"air_pressure", "positive"=>"down"))
        v_lev[:] = target_levels

    elseif file_type == :raw && vert_coord !== nothing
        # Isca-style native vertical metadata. These one-dimensional pressure
        # axes are approximate reference values; exact column pressure is
        # reconstructed from interface pk, bk, and the output-record ps.
        defDim(ds, "pfull", mesh.nd)
        defDim(ds, "phalf", mesh.nd + 1)

        reference_ps = fill(Float64(vert_coord.p_ref), 1, 1)
        reference_pfull = zeros(Float64, 1, 1, mesh.nd)
        Compute_Pressure_Grid!(reference_pfull, vert_coord.ak, vert_coord.bk, reference_ps)
        reference_phalf = vert_coord.ak .+ vert_coord.bk .* vert_coord.p_ref

        v_pfull = defVar(ds, "pfull", Float64, ("pfull",), attrib=Dict(
            "units"=>"hPa", "long_name"=>"approx full pressure level",
            "positive"=>"down", "axis"=>"Z"))
        v_pfull[:] = vec(reference_pfull) .* 0.01

        v_phalf = defVar(ds, "phalf", Float64, ("phalf",), attrib=Dict(
            "units"=>"hPa", "long_name"=>"approx half pressure level",
            "positive"=>"down", "axis"=>"Z"))
        v_phalf[:] = reference_phalf .* 0.01

        defVar(ds, "pk", Float64, ("phalf",), attrib=Dict(
            "units"=>"Pa", "long_name"=>"vertical coordinate pressure values"))[:] = vert_coord.ak
        defVar(ds, "bk", Float64, ("phalf",), attrib=Dict(
            "units"=>"1", "long_name"=>"vertical coordinate sigma values"))[:] = vert_coord.bk
    end
    # --- Coordinate Variables & Vertical Metadata --- #

    # --- Data Variables --- #
    for sym in requested_vars
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym]

        dims = if meta.dims == 3
            (file_type == :plev) ? ("lon", "lat", "plev", "time") : ("lon", "lat", "pfull", "time")
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
    start_time, ::Any, requested_vars;
    filename::String = "output.nc",
    do_plev_output::Bool = false,
    pressure_levels::Vector{Float64} = Float64[],
    day_to_sec::Int64 = 86400,
    output_interval::Int64 = 86400,
    spinup_day::Float64 = 0.0,
    model_mode::Symbol = :PrimitiveEquation,
    institute::String = "My Research Lab",
    experiment_id::String = "Simulation_v1"
)
    day_to_sec > 0 || throw(ArgumentError("day_to_sec must be positive"))
    output_interval > 0 || throw(ArgumentError("output_interval must be positive"))
    spinup_day >= 0 || throw(ArgumentError("spinup_day must be non-negative"))
    if do_plev_output
        !isempty(pressure_levels) || throw(
            ArgumentError("pressure-level output requires at least one pressure level"),
        )
        all(p -> isfinite(p) && p > 0, pressure_levels) || throw(
            ArgumentError("pressure levels must be positive and finite"),
        )
    end

    # --- Instantiate the Mode Type --- #
    mode_obj    = Mode_Factory(model_mode)
    if do_plev_output && !isa(mode_obj, PrimitiveEquationMode)
        throw(ArgumentError("pressure-level output is only available for PrimitiveEquation"))
    end
    nλ, nθ, nd  = mesh.nλ, mesh.nθ, mesh.nd
    spinup_time = start_time + Int64(spinup_day * day_to_sec)
    active_symbols = unique(copy(requested_vars))

    # Derive base filename (strip ".nc" extension)
    base_fn = splitext(filename)[1]

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
    vert_ak    = Float64[]
    vert_bk    = Float64[]
    # --- Instantiate the Mode Type --- #

    # --- 3D Specific Initialization --- #
    if isa(mode_obj, PrimitiveEquationMode)
        # Validation
        if !(:ps in active_symbols); push!(active_symbols, :ps); end
        if (:z in active_symbols) && !(:t in active_symbols); error("Need :t for :z"); end

        if vert_coord === nothing
             error("Output_Manager: vert_coord is required for :PrimitiveEquation mode")
        end

        # Allocation (Only happens for 3D)
        n_plev     = length(pressure_levels)
        p3d_buf    = zeros(Float64, nλ, nθ, nd)
        interp_buf = zeros(Float64, nλ, nθ, max(n_plev, 1))

        vert_ak = copy(vert_coord.ak)
        vert_bk = copy(vert_coord.bk)
    end

    # Fetch variable info
    var_info_map = Get_Var_Info(Val(mode_symbol(mode_obj)))
    unknown_symbols = setdiff(active_symbols, keys(var_info_map))
    isempty(unknown_symbols) || throw(
        ArgumentError("unsupported output variables for $(mode_symbol(mode_obj)): $unknown_symbols"),
    )
    # --- 3D Specific Initialization --- #

    # --- Initialize Accumulators --- #
    # acc_raw: always allocated (native grid accumulator for all model types)
    acc_raw = Dict{Symbol, Array{Float64}}()
    for sym in active_symbols
        haskey(var_info_map, sym) || continue
        meta = var_info_map[sym]
        acc_raw[sym] = (meta.dims == 3) ? zeros(Float64, nλ, nθ, nd) : zeros(Float64, nλ, nθ)
    end

    # acc_plev: only allocated when pressure-level output is requested (3D only)
    acc_plev = (do_plev_output && isa(mode_obj, PrimitiveEquationMode)) ?
               deepcopy(acc_raw) : Dict{Symbol, Array{Float64}}()
    # --- Initialize Accumulators --- #

    # --- Initialize Files --- #
    # ds_raw: always created for all model types (primary/native output)
    raw_file_type = isa(mode_obj, PrimitiveEquationMode) ? :raw : :simple_2d
    ds_raw = _Init_Single_File(
        "$(base_fn)_t$(start_time).nc", mesh,
        var_info_map, active_symbols, raw_file_type;
        vert_coord=vert_coord, global_meta=global_meta
    )

    # ds_plev: only for 3D with do_plev_output = true
    ds_plev = nothing
    if do_plev_output && isa(mode_obj, PrimitiveEquationMode)
        ds_plev = _Init_Single_File(
            "$(base_fn)_t$(start_time)_plev.nc", mesh,
            var_info_map, active_symbols, :plev, pressure_levels;
            vert_coord=vert_coord, global_meta=global_meta
        )
    end
    # --- Initialize Files --- #

    return Output_Manager(
        nλ, nθ, nd, atmo_data, mode_obj,
        do_plev_output, pressure_levels, log.(max.(pressure_levels, 1.0)),
        p3d_buf, interp_buf, vert_ak, vert_bk,
        day_to_sec, start_time, start_time, spinup_time, output_interval, 0, 1,
        ds_plev, ds_raw, acc_plev, acc_raw, active_symbols, var_info_map,
        base_fn, collect(mesh.λc), collect(mesh.θc), vert_coord, global_meta
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

        # acc_raw: always accumulate (native sigma/hybrid grid → ds_raw)
        if haskey(mgr.acc_raw, sym)
            _accumulate_buffer!(mgr.acc_raw[sym], raw_data)
        end

        # acc_plev: only accumulate when pressure-level output is requested (→ ds_plev)
        if mgr.do_plev_output && haskey(mgr.acc_plev, sym)
            _accumulate_buffer!(mgr.acc_plev[sym], raw_data)
        end
    end
end

# --- Specialized Accumulation: 2D Modes ---
function _accumulate_core!(::AbstractModelMode, mgr, live_data)
    # Generic fallback for Barotropic/ShallowWater: use acc_raw
    for sym in mgr.active_symbols
        if !haskey(live_data, sym); continue; end

        if haskey(mgr.acc_raw, sym)
            mgr.acc_raw[sym] .+= live_data[sym]
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

    t_val = manager.current_time / manager.day_to_sec

    # Write time index to whichever files are open
    if manager.ds_raw  !== nothing
        manager.ds_raw["time"][manager.output_index]  = t_val
    end
    if manager.ds_plev !== nothing
        manager.ds_plev["time"][manager.output_index] = t_val
    end

    # Dispatch to specific writing logic
    _write_core!(manager.mode, manager)

    manager.sample_counter = 0
    manager.output_index += 1
end

# --- Writer: Primitive Equation ---
function _write_core!(::PrimitiveEquationMode, mgr)
    t_idx        = mgr.output_index
    sample_count = Float64(mgr.sample_counter)
    var_info_map = mgr.var_info_map

    # =========================================================================
    # PRIMARY: Native sigma/hybrid output (acc_raw → ds_raw, always)
    # =========================================================================
    for sym in keys(mgr.acc_raw)
        haskey(var_info_map, sym) || continue
        meta     = var_info_map[sym]
        nc_var   = mgr.ds_raw[meta.nc_name]
        raw_mean = mgr.acc_raw[sym]
        raw_mean ./= sample_count

        if ndims(nc_var) == 4  # (lon, lat, pfull, time)
            nc_var[:, :, :, t_idx] = raw_mean
        elseif ndims(nc_var) == 3  # (lon, lat, time)
            nc_var[:, :, t_idx] = (ndims(raw_mean) == 3) ? view(raw_mean, :, :, 1) : raw_mean
        end
        fill!(mgr.acc_raw[sym], 0.0)
    end
    NCDatasets.sync(mgr.ds_raw)

    # =========================================================================
    # OPTIONAL: Pressure-level interpolated output (acc_plev → ds_plev)
    # =========================================================================
    if mgr.do_plev_output && mgr.ds_plev !== nothing
        for accumulator in values(mgr.acc_plev)
            accumulator ./= sample_count
        end
        ps_avg_2d = mgr.acc_plev[:ps]
        Compute_Pressure_Grid!(mgr.p3d_buffer, mgr.vert_ak, mgr.vert_bk, ps_avg_2d)

        for sym in keys(mgr.acc_plev)
            haskey(var_info_map, sym) || continue
            meta        = var_info_map[sym]
            nc_var      = mgr.ds_plev[meta.nc_name]
            native_mean = mgr.acc_plev[sym]

            if ndims(nc_var) == 4  # (lon, lat, plev, time)
                t_ref = get(mgr.acc_plev, :t, nothing)
                Interpolate_Field!(mgr.interp_buffer, native_mean, mgr.p3d_buffer, ps_avg_2d,
                                   mgr.log_targets, sym, mgr.atmo_data, t_ref)
                nc_var[:, :, :, t_idx] = mgr.interp_buffer
            else
                nc_var[:, :, t_idx] = (ndims(native_mean) == 3) ? view(native_mean, :, :, 1) : native_mean
            end
        end
        for accumulator in values(mgr.acc_plev)
            fill!(accumulator, 0.0)
        end
        NCDatasets.sync(mgr.ds_plev)
    end
end

# --- Writer: 2D Modes (Barotropic/ShallowWater) ---
function _write_core!(::AbstractModelMode, mgr)
    t_idx        = mgr.output_index
    sample_count = Float64(mgr.sample_counter)
    var_info_map = mgr.var_info_map

    for sym in keys(mgr.acc_raw)
        if !haskey(var_info_map, sym); continue; end
        nc_name = var_info_map[sym].nc_name
        accumulator = mgr.acc_raw[sym]
        accumulator ./= sample_count
        mgr.ds_raw[nc_name][:, :, t_idx] = accumulator
        fill!(accumulator, 0.0)
    end
    NCDatasets.sync(mgr.ds_raw)
end

# ==============================================================================
# 7. Chunk Rotation
# ==============================================================================

"""
    _open_chunk_files!(manager, chunk_time)

Internal helper: open a new pair of NC output files for the given chunk start time.
Uses a NamedTuple mesh proxy so the full mesh object need not be stored in the manager.
"""
function _open_chunk_files!(manager::Output_Manager, chunk_time::Int64)
    var_info_map = manager.var_info_map

    # NamedTuple acts as a mesh proxy (supports dot-notation like a struct)
    mesh_proxy = (
        nλ = manager.nλ,
        nθ = manager.nθ,
        nd = manager.nd,
        λc = manager.mesh_λc,
        θc = manager.mesh_θc
    )

    # ds_raw: always (primary output for all model types)
    raw_file_type = isa(manager.mode, PrimitiveEquationMode) ? :raw : :simple_2d
    manager.ds_raw = _Init_Single_File(
        "$(manager.base_filename)_t$(chunk_time).nc", mesh_proxy,
        var_info_map, manager.active_symbols, raw_file_type;
        vert_coord=manager.vert_coord_cache, global_meta=manager.global_meta
    )

    # ds_plev: optional (3D with do_plev_output = true only)
    if manager.do_plev_output && isa(manager.mode, PrimitiveEquationMode)
        manager.ds_plev = _Init_Single_File(
            "$(manager.base_filename)_t$(chunk_time)_plev.nc", mesh_proxy,
            var_info_map, manager.active_symbols, :plev, manager.target_levels;
            vert_coord=manager.vert_coord_cache, global_meta=manager.global_meta
        )
    end
end

"""
    Rotate_NC_Chunk!(manager, new_chunk_time)

Flush any remaining data, close the current NC output files, and open fresh
time-stamped files for the next chunk. Called at each saving_frequency boundary
in the driver, coordinated with JLD2 checkpoint writes.
"""
function Rotate_NC_Chunk!(manager::Output_Manager, new_chunk_time::Int64)
    # Flush any partially-accumulated data into the closing chunk
    if manager.sample_counter > 0
        Flush_to_Disk!(manager)
    end

    # Close current chunk files
    if manager.ds_plev !== nothing
        close(manager.ds_plev)
        manager.ds_plev = nothing
    end
    if manager.ds_raw !== nothing
        close(manager.ds_raw)
        manager.ds_raw = nothing
    end

    # Open next chunk files
    _open_chunk_files!(manager, new_chunk_time)

    # Reset time index for the new file
    manager.output_index = 1

    @info "NC output rotated: new chunk at t=$(new_chunk_time)s → $(manager.base_filename)_t$(new_chunk_time).nc"
end

# ==============================================================================
# 8. Finalization
# ==============================================================================

function Finalize_Output!(manager::Output_Manager)
    if manager.sample_counter > 0
        Flush_to_Disk!(manager)
    end
    if manager.ds_plev !== nothing; close(manager.ds_plev); end
    if manager.ds_raw  !== nothing; close(manager.ds_raw);  end
end

end
