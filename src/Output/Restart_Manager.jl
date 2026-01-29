module Restart_Manager_Module

using JLD2
using ..Dyn_Data_Module

export Restart_Manager, Write_Restart_File, Load_Restart_File!


struct Restart_Manager
    output_dir::String
    restart_frequency::Int64 
end



"""
    Write_Restart_File(manager, dyn_data, current_time)
"""
function Write_Restart_File(manager::Restart_Manager, dyn_data::Dyn_Data, current_time::Int64)
    # Ensure directory exists right before writing (lazy creation)
    if !isdir(manager.output_dir)
        mkpath(manager.output_dir)
    end
    
    filename = joinpath(manager.output_dir, "restart_t$(current_time).jld2")
    temp_filename = filename * ".tmp"

    jldsave(temp_filename; 
        dyn_data_state = dyn_data, 
        saved_time = current_time
    )

    mv(temp_filename, filename; force=true)
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

end