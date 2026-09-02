using NCDatasets
using Printf
using JGCM

import JGCM.Atmo_Data_Module: Atmo_Data
import JGCM.Output_Mappings_Module: Get_Var_Info
import JGCM.Vertical_Interpolation_Module:
    Pressure_Interpolation_Cache, Prepare_Interpolation!, Apply_Interpolation!

"""
    Interpolate_File(input_path, output_path, target_levels_Pa; var_names)
"""
function Interpolate_File(
    input_path::String,
    output_path::String,
    target_levels::Vector{Float64};
    var_names::Vector{Symbol} = [:u, :v, :t, :z, :q, :vor, :div],
)
    isempty(target_levels) && throw(ArgumentError("at least one pressure level is required"))
    all(p -> isfinite(p) && p > 0.0, target_levels) ||
        throw(ArgumentError("pressure levels must be positive and finite"))

    var_info = Get_Var_Info(Val(:PrimitiveEquation))
    requested_symbols = unique(var_names)
    unknown_symbols = setdiff(requested_symbols, keys(var_info))
    isempty(unknown_symbols) || throw(
        ArgumentError("unsupported primitive-equation variables: $unknown_symbols"),
    )
    requested_variables = [(sym, var_info[sym].nc_name) for sym in requested_symbols]

    @printf("--- Starting Post-Processing ---\n")
    @printf("Input:  %s\n", input_path)
    @printf("Output: %s\n", output_path)
    @printf("Levels: %s Pa\n", target_levels)

    # 1. Open Input File
    ds_in = NCDataset(input_path, "r")

    # 2. Read Isca-style interface coefficients. Arithmetic midpoint
    # coefficients are intentionally unsupported because they do not reproduce
    # the model's Simmons--Burridge full-level pressure.
    haskey(ds_in, "pk") && haskey(ds_in, "bk") ||
        error("Vertical grid coefficients not found. Need interface variables 'pk' and 'bk'.")
    haskey(ds_in.dim, "pfull") ||
        error("Full-level dimension 'pfull' is missing from the input file.")
    ak = Float64.(ds_in["pk"][:])
    bk = Float64.(ds_in["bk"][:])
    @printf("Grid: Detected Isca-style interface coefficients (pk/bk)\n")

    # 3. Setup Physics Context. Native files record the constants used by the
    # model so offline geopotential-height interpolation matches online output.
    gravity = haskey(ds_in.attrib, "gravity") ? Float64(ds_in.attrib["gravity"]) : 9.80
    rdgas = haskey(ds_in.attrib, "dry_air_gas_constant") ?
            Float64(ds_in.attrib["dry_air_gas_constant"]) : 287.04
    isfinite(gravity) && gravity > 0.0 ||
        throw(DomainError(gravity, "gravity metadata must be finite and positive"))
    isfinite(rdgas) && rdgas > 0.0 || throw(
        DomainError(rdgas, "dry-air gas-constant metadata must be finite and positive"),
    )
    phys = Atmo_Data(
        "PostProcess",
        1,
        1,
        1,
        false,
        false,
        false,
        false,
        [0.0];
        radius = 6371.0e3,
        grav = gravity,
        rdgas = rdgas,
    )

    # 4. Prepare Output File
    if isfile(output_path)
        rm(output_path)
    end
    ds_out = NCDataset(output_path, "c")

    # Define Dimensions
    defDim(ds_out, "lon", ds_in.dim["lon"])
    defDim(ds_out, "lat", ds_in.dim["lat"])
    defDim(ds_out, "time", Inf)
    defDim(ds_out, "plev", length(target_levels))

    # Copy Coordinates (lon, lat, time)
    for coord in ["lon", "lat", "time"]
        if haskey(ds_in, coord)
            v_dst = defVar(ds_out, coord, Float64, (coord,), attrib = ds_in[coord].attrib)
            coord_values = ds_in[coord][:]
            v_dst[1:length(coord_values)] = coord_values
        end
    end

    # Define Pressure Levels variable
    v_plev = defVar(
        ds_out,
        "plev",
        Float64,
        ("plev",),
        attrib = Dict("units" => "Pa", "axis" => "Z"),
    )
    v_plev[:] = target_levels

    # [NEW] Define and Copy Surface Pressure (PS)
    # We always want PS in the output for reference, but keeping its original 2D structure.
    if haskey(ds_in, "ps")
        defVar(ds_out, "ps", Float64, ("lon", "lat", "time"), attrib = ds_in["ps"].attrib)
    end

    # Define 3D Interpolated Variables
    for (_, nc_name) in requested_variables
        if haskey(ds_in, nc_name) && ndims(ds_in[nc_name]) == 4
            defVar(
                ds_out,
                nc_name,
                Float64,
                ("lon", "lat", "plev", "time"),
                attrib = ds_in[nc_name].attrib,
            )
        end
    end

    # 5. Initialize Buffers
    nλ = ds_in.dim["lon"]
    nθ = ds_in.dim["lat"]
    nd = ds_in.dim["pfull"]
    length(ak) == nd + 1 && length(bk) == nd + 1 || throw(
        DimensionMismatch("pk and bk must each contain pfull + 1 interface values"),
    )
    n_tgt = length(target_levels)
    log_targets = log.(target_levels)

    interpolation_cache = Pressure_Interpolation_Cache(nλ, nθ, nd, n_tgt)
    interp_buffer = zeros(Float64, nλ, nθ, n_tgt)

    temperature_name = var_info[:t].nc_name
    needs_temperature = any(
        var_sym == :z && haskey(ds_in, nc_name) && ndims(ds_in[nc_name]) == 4
        for (var_sym, nc_name) in requested_variables
    )

    # 6. Time Loop
    n_times = ds_in.dim["time"]

    for t = 1:n_times
        @printf("Processing Time Step: %d / %d\r", t, n_times)

        # A. Read & Copy Surface Pressure
        if !haskey(ds_in, "ps")
            error("Variable 'ps' is missing in input file.")
        end
        ps_slice = ds_in["ps"][:, :, t]

        # Write ps to output (2D copy)
        ds_out["ps"][:, :, t] = ps_slice

        # B. Prepare pressure geometry once and reuse it for every 3D field,
        # following Isca's pres_interp_type setup/apply split.
        Prepare_Interpolation!(interpolation_cache, ak, bk, ps_slice, log_targets)

        # C. Read Temperature (for Hydrostatic calc if needed)
        t_ref = nothing
        if needs_temperature && haskey(ds_in, temperature_name)
            t_ref = ds_in[temperature_name][:, :, :, t]
        end

        # D. Process 3D Variables
        for (var_sym, nc_name) in requested_variables
            if !haskey(ds_in, nc_name)
                continue
            end

            if ndims(ds_in[nc_name]) != 4
                continue
            end

            # If height needs temperature, reuse the already-loaded slice when
            # temperature itself is also requested.
            raw_data = var_sym == :t && !isnothing(t_ref) ?
                       t_ref : ds_in[nc_name][:, :, :, t]

            Apply_Interpolation!(
                interp_buffer,
                raw_data,
                interpolation_cache,
                var_sym,
                phys,
                t_ref,
            )

            ds_out[nc_name][:, :, :, t] = interp_buffer
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
        println(
            "Usage: julia --project post_processing/Interpolator.jl <input.nc> <output.nc> [levels...]",
        )
        exit(1)
    end

    input_fn = ARGS[1]
    output_fn = ARGS[2]

    if length(ARGS) > 2
        p_levels = parse.(Float64, ARGS[3:end])
    else
        println("No levels provided. Defaulting to [85000, 50000, 20000]")
        p_levels = [85000.0, 50000.0, 20000.0]
    end

    Interpolate_File(input_fn, output_fn, p_levels)
end
