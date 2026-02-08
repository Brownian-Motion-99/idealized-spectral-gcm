using NCDatasets
using Base.Threads
using Printf
using JGCM 

import JGCM.Atmo_Data_Module: Atmo_Data
import JGCM.Vertical_Interpolation_Module: Interpolate_Field!

"""
    Interpolate_File(input_path, output_path, target_levels_Pa; var_names)
"""
function Interpolate_File(
    input_path::String, 
    output_path::String, 
    target_levels::Vector{Float64};
    var_names::Vector{Symbol}=[:u, :v, :t, :z, :q, :vor, :div]
)
    @printf("--- Starting Post-Processing ---\n")
    @printf("Input:  %s\n", input_path)
    @printf("Output: %s\n", output_path)
    @printf("Levels: %s Pa\n", target_levels)
    
    # 1. Open Input File
    ds_in = NCDataset(input_path, "r")
    
    # 2. Reconstruct Vertical Grid
    ak_mid, bk_mid = Float64[], Float64[]

    if haskey(ds_in, "hyam") && haskey(ds_in, "hybm")
        @printf("Grid: Detected Midpoint Coefficients (hyam/hybm)\n")
        ak_mid = ds_in["hyam"][:]
        bk_mid = ds_in["hybm"][:]
        
    elseif haskey(ds_in, "hyai") && haskey(ds_in, "hybi")
        @printf("Grid: Detected Interface Coefficients (hyai/hybi) -> Calculating Midpoints\n")
        ak_int = ds_in["hyai"][:]
        bk_int = ds_in["hybi"][:]
        ak_mid = 0.5 .* (ak_int[1:end-1] .+ ak_int[2:end])
        bk_mid = 0.5 .* (bk_int[1:end-1] .+ bk_int[2:end])
    else
        error("Vertical grid coefficients not found. Need ('hyam', 'hybm') or ('hyai', 'hybi').")
    end
    
    # 3. Setup Physics Context
    phys = Atmo_Data(
        "PostProcess", 
        1, 1, 1, 
        false, false, false, false, 
        [0.0]; 
        radius = 6371.0e3
    )

    # 4. Prepare Output File
    if isfile(output_path); rm(output_path); end
    ds_out = NCDataset(output_path, "c")
    
    # Define Dimensions
    defDim(ds_out, "lon",  ds_in.dim["lon"])
    defDim(ds_out, "lat",  ds_in.dim["lat"])
    defDim(ds_out, "time", Inf) 
    defDim(ds_out, "plev", length(target_levels))
    
    # Copy Coordinates (lon, lat, time)
    for coord in ["lon", "lat", "time"]
        if haskey(ds_in, coord)
            v_dst = defVar(ds_out, coord, Float64, (coord,), attrib=ds_in[coord].attrib)
            v_dst[:] = ds_in[coord][:]
        end
    end
    
    # Define Pressure Levels variable
    v_plev = defVar(ds_out, "plev", Float64, ("plev",), attrib=Dict("units"=>"Pa", "axis"=>"Z"))
    v_plev[:] = target_levels

    # [NEW] Define and Copy Surface Pressure (PS)
    # We always want PS in the output for reference, but keeping its original 2D structure.
    if haskey(ds_in, "ps")
        defVar(ds_out, "ps", Float64, ("lon", "lat", "time"), attrib=ds_in["ps"].attrib)
    end

    # Define 3D Interpolated Variables
    for var_sym in var_names
        var_str = String(var_sym)
        if haskey(ds_in, var_str) && ndims(ds_in[var_str]) == 4
            defVar(ds_out, var_str, Float64, ("lon", "lat", "plev", "time"), 
                   attrib=ds_in[var_str].attrib)
        end
    end

    # 5. Initialize Buffers
    nλ = ds_in.dim["lon"]
    nθ = ds_in.dim["lat"]
    nd = ds_in.dim["lev"]
    n_tgt = length(target_levels)
    log_targets = log.(target_levels)
    
    p3d_buffer    = zeros(Float64, nλ, nθ, nd)
    interp_buffer = zeros(Float64, nλ, nθ, n_tgt)
    
    # 6. Time Loop
    n_times = ds_in.dim["time"]

    for t in 1:n_times
        @printf("Processing Time Step: %d / %d\r", t, n_times)
        
        # A. Read & Copy Surface Pressure
        if !haskey(ds_in, "ps")
            error("Variable 'ps' is missing in input file.")
        end
        ps_slice = ds_in["ps"][:, :, t]
        
        # Write ps to output (2D copy)
        ds_out["ps"][:, :, t] = ps_slice
        
        # B. Reconstruct 3D Pressure Field
        @threads for j in 1:nθ
            for i in 1:nλ
                ps_val = ps_slice[i, j]
                for k in 1:nd
                    p3d_buffer[i, j, k] = ak_mid[k] + bk_mid[k] * ps_val
                end
            end
        end
        
        # C. Read Temperature (for Hydrostatic calc if needed)
        t_ref = nothing
        if haskey(ds_in, "t")
            t_ref = ds_in["t"][:, :, :, t]
        end

        # D. Process 3D Variables
        for var_sym in var_names
            var_str = String(var_sym)
            if !haskey(ds_in, var_str); continue; end
            
            if ndims(ds_in[var_str]) != 4
                continue 
            end
            
            raw_data = ds_in[var_str][:, :, :, t]
            
            Interpolate_Field!(
                interp_buffer,
                raw_data,
                p3d_buffer,
                ps_slice,    
                log_targets,
                var_sym,
                phys,
                t_ref        
            )
            
            ds_out[var_str][:, :, :, t] = interp_buffer
        end
        
        if t % 10 == 0
            NCDatasets.sync(ds_out)
        end
    end
    
    close(ds_in)
    close(ds_out)
    println("\nDone! Output saved to: $output_path")
end

# CLI Entry Point
# Usage: julia --project= post_processing/Interpolator.jl $INPUT_FILE $OUTPUT_FILE [levels...]
if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 2
        println("Usage: julia --project post_processing/Interpolator.jl <input.nc> <output.nc> [levels...]")
        exit(1)
    end

    input_fn  = ARGS[1]
    output_fn = ARGS[2]
    
    if length(ARGS) > 2
        p_levels = parse.(Float64, ARGS[3:end])
    else
        println("No levels provided. Defaulting to [85000, 50000, 20000]")
        p_levels = [85000.0, 50000.0, 20000.0]
    end

    Interpolate_File(input_fn, output_fn, p_levels) 
end