# Export the binary staging inputs from Build_Experiment_LRF.py and validate
# every result with the actual model loader/kernel and current T42L20 grid.
# julia --project=. --compiled-modules=no --threads=4 \
#   post_processing/Write_Experiment_LRF.jl STAGING_DIRECTORY
using JGCM, JLD2, TOML, SHA, Printf

function read_f64(directory, name, dimensions)
    path = joinpath(directory, name * ".f64")
    values = Vector{Float64}(undef, prod(dimensions))
    open(path, "r") do io
        read!(io, values)
        eof(io) || error("Unexpected trailing bytes in $path")
    end
    return reshape(values, dimensions)
end

function main()
    length(ARGS) == 1 || error("Usage: Write_Experiment_LRF.jl STAGING_DIRECTORY")
    root = abspath(ARGS[1])
    mesh = Spectral_Spherical_Mesh(42, 43, 128, 64, 20, 6.371e6)
    latitude_model = rad2deg.(mesh.θc)
    vert = Vert_Coordinate(128, 64, 20, "even_sigma",
        "simmons_and_burridge", "second_centered_wts")
    nominal_pfull = zeros(1, 1, 20)
    JGCM.Vertical_Interpolation_Module.Compute_Pressure_Grid!(
        nominal_pfull, vert.ak, vert.bk, fill(vert.p_ref, 1, 1))
    nominal_phalf = (vert.ak .+ vert.bk .* vert.p_ref) ./ 100
    results = Dict{String,Any}(
        "julia_version" => string(VERSION), "threads" => Threads.nthreads(),
        "writer_sha256" => bytes2hex(sha256(read(@__FILE__))),
        "lrf_kernel_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "../src/Physics/LRF.jl")))),
    )
    for experiment in ("ctrl_BM", "sst1.0_BM", "sst2.5_BM")
        directory = joinpath(root, "export_inputs", experiment)
        metadata = TOML.parsefile(joinpath(directory, "metadata.toml"))
        nd, nlat = metadata["nlevel"], metadata["nlatitude"]
        (nd, nlat) == (20, 64) || error("Expected T42L20 dimensions")
        B_by_lat = read_f64(directory, "B_by_lat", (nd, nd, nlat))
        chi_reference = read_f64(directory, "chi_reference", (nd, nlat))
        latitude = vec(read_f64(directory, "latitude", (nlat,)))
        pfull = vec(read_f64(directory, "pfull", (nd,)))
        phalf = vec(read_f64(directory, "phalf", (nd + 1,)))
        validation_q = read_f64(directory, "validation_q", (1, nlat, nd))
        validation_tendency_K_s = read_f64(directory, "validation_tendency_K_s", (1, nlat, nd))
        maximum(abs, latitude .- latitude_model) <= 1e-10 || error("Wrong model latitudes")
        pfull_error = maximum(abs, pfull .- vec(nominal_pfull) ./ 100)
        phalf_error = maximum(abs, phalf .- nominal_phalf)
        pfull_error <= 1e-9 && phalf_error <= 1e-9 || error("Wrong nominal model pressures")
        hdf_path = joinpath(root, experiment * ".h5")
        source_sha256 = bytes2hex(sha256(read(hdf_path)))
        source_sha256 == metadata["source_sha256"] || error("HDF5 source hash mismatch")
        output = joinpath(root, experiment * ".jld2")
        ispath(output) && error("Refusing to overwrite $output")
        JLD2.jldsave(output;
            scheme = "regularized_log_latitude_v1",
            formulation = "Q_K_day[:,j] = alpha * B_by_lat[:,:,j] * (log(q[:,j] + q0) - chi_reference[:,j])",
            B_by_lat, chi_reference,
            q0_kg_kg = metadata["q0_kg_kg"], alpha = metadata["alpha"],
            latitude, pfull, phalf, source_sha256,
            source_file = experiment * ".h5", experiment,
            metadata_json = metadata["metadata_json"],
            vertical_order = "top_to_bottom",
            B_dimensions = "response_level,perturbation_level,latitude",
            B_units = "K day-1", humidity_units = "kg kg-1",
            validation_q, validation_tendency_K_s,
        )

        # Full model grid compatibility, and the embedded independent Python
        # calculation on one longitude at every model latitude.
        full_state = Load_LRF_State(output, 128, 64, 20; latitude = latitude_model)
        full_state isa Latitude_LRF_State || error("Wrong loaded scheme")
        state = Load_LRF_State(output, 1, 64, 20; latitude = latitude_model)
        actual = similar(validation_q)
        before = copy(validation_q)
        LRF!(state, validation_q, actual, 86400)
        validation_q == before || error("LRF mutated humidity")
        cross_language_error = maximum(abs, actual .- validation_tendency_K_s)
        scale = maximum(abs, validation_tendency_K_s)
        cross_language_error <= 5e-14 * max(scale, 1e-12) || error("Python/Julia heating mismatch")

        reference = exp.(chi_reference) .- state.q0
        minimum(reference) >= 0 || error("Nonphysical humidity reference")
        reference_grid = permutedims(reshape(reference, nd, nlat, 1), (3, 2, 1))
        LRF!(state, reference_grid, actual, 86400)
        reference_error = maximum(abs, actual)
        reference_error <= 1e-17 || error("Nonzero reference heating")
        LRF!(state, zeros(1, nlat, nd), actual, 86400)
        all(isfinite, actual) || error("Nonfinite zero-humidity response")

        results[experiment] = Dict(
            "file" => experiment * ".jld2", "scheme" => "regularized_log_latitude_v1",
            "sha256" => bytes2hex(sha256(read(output))), "source_sha256" => source_sha256,
            "cross_language_max_error_K_s" => cross_language_error,
            "maximum_validation_tendency_K_s" => scale,
            "reference_max_error_K_s" => reference_error,
            "nominal_pfull_max_error_hPa" => pfull_error,
            "nominal_phalf_max_error_hPa" => phalf_error,
            "q0_kg_kg" => state.q0, "alpha" => state.alpha,
            "grid_dimensions" => [128, 64, 20], "validated" => true,
        )
        @printf("Validated %s: Python/Julia error %.3e K/s, reference %.3e K/s\n",
            output, cross_language_error, reference_error)
        flush(stdout)
    end
    open(joinpath(root, "julia_validation.toml"), "w") do io
        TOML.print(io, results)
    end
end

main()
