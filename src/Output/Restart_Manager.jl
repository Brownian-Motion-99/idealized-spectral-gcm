module Restart_Manager_Module

using JLD2
using ..Dyn_Data_Module

export Restart_Manager, Write_Restart_File, Load_Restart_File!, Cleanup_Old_Restarts


struct Restart_Manager
    output_dir::String
    restart_frequency::Int64

    function Restart_Manager(output_dir::String, frequency::Int64)
        if frequency > 0 && !isdir(output_dir)
            mkpath(output_dir)
        end
        return new(output_dir, frequency)
    end
end



"""
    Write_Restart_File(manager, dyn_data, current_time)
"""
function Write_Restart_File(
    manager::Restart_Manager,
    dyn_data::Dyn_Data,
    current_time::Int64,
)
    # Ensure directory exists right before writing (lazy creation)
    if !isdir(manager.output_dir)
        mkpath(manager.output_dir)
    end

    filename = joinpath(manager.output_dir, "restart_t$(current_time).jld2")
    temp_filename = filename * ".tmp"

    jldsave(temp_filename; dyn_data_state = dyn_data, saved_time = current_time)

    mv(temp_filename, filename; force = true)
    @info "Checkpoint saved: $filename"
end



"""
    Load_Restart_File!(dyn_data, filename)
"""
function Load_Restart_File!(dyn_data::Dyn_Data, filename::String)
    if !isfile(filename)
        error("Restart file not found: $filename")
    end

    @info "Loading warm start from: $filename"

    loaded_file = load(filename)
    loaded_struct = loaded_file["dyn_data_state"]
    saved_time = loaded_file["saved_time"]

    for name in fieldnames(Dyn_Data)
        src_field = getfield(loaded_struct, name)
        dest_field = getfield(dyn_data, name)

        if isa(src_field, Array)
            copyto!(dest_field, src_field)
        else
            setfield!(dyn_data, name, src_field)
        end
    end

    return saved_time
end



"""
    Cleanup_Old_Restarts(manager, keep_last_n, protect_file)
    
Deletes older restart files, BUT ensures the file used to start the run 
(protect_file) is never deleted.
"""
function Cleanup_Old_Restarts(
    manager::Restart_Manager,
    keep_last_n::Int = 3,
    protect_file::String = "",
)
    files = readdir(manager.output_dir, join = true)
    restart_files = filter(f -> occursin("restart_t", f) && endswith(f, ".jld2"), files)

    # Helper to extract time
    function get_time_step(f)
        m = match(r"restart_t(\d+).jld2", f)
        return isnothing(m) ? 0 : parse(Int64, m[1])
    end

    sort!(restart_files, by = get_time_step)

    if length(restart_files) > keep_last_n
        files_to_delete = restart_files[1:end-keep_last_n]

        for f in files_to_delete
            # SAFETY CHECK: Do not delete if it matches the start file
            if !isempty(protect_file) && abspath(f) == abspath(protect_file)
                @info "Skipping cleanup for protected start file: $f"
                continue
            end

            rm(f; force = true)
            @info "Disk cleanup: Deleted old checkpoint $f"
        end
    end
end

end
