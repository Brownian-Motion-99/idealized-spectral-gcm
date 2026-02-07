module Experiment_Configuration

import Base

export Model_Config

"""
    Model_Config

The blueprint for a JGCM simulation. 
"""
Base.@kwdef struct Model_Config
    # -------------------------------------------------------
    # 1. Experiment Identity
    # -------------------------------------------------------
    name::String
    model_type::Symbol  # :Barotropic, :Shallow_Water, :PrimitiveEquation

    # -------------------------------------------------------
    # 2. Resolution & Geometry
    # -------------------------------------------------------
    num_fourier::Int    
    nθ::Int             
    nd::Int             
    
    radius::Float64     
    omega::Float64
    grav::Float64

    # -------------------------------------------------------
    # 3. Vertical Coordinate
    # -------------------------------------------------------
    # Options: "even_sigma", "uneven_sigma", "hybrid", "simmons_and_burridge", "mcm", "v197"
    vert_coord_option::Any    
    # Options: "simmons_and_burridge".
    vert_difference_option::Any 
    # Options: "second_centered_wts", "second_centered".
    vert_ref_level_option::Any  

    # -------------------------------------------------------
    # 4. Time Integration & Planet Settings
    # -------------------------------------------------------
    Δt::Int64
    end_time::Int64
    spinup_day::Float64 = 0.0
    day_to_sec::Int64
    
    damping_order::Int      
    damping_coef::Float64   
    robert_coef::Float64
    implicit_coef::Float64  

    # -------------------------------------------------------
    # 5. Restart Configuration
    # -------------------------------------------------------
    is_restart::Bool         = false       # Toggle warm start
    restart_file::String     = ""          # Path to the .jld2 file to load
    restart_frequency::Int64 = 86400       # How often to save checkpoints (seconds)

    # -------------------------------------------------------
    # 6. Composition & Physics
    # -------------------------------------------------------   
    moisture_processes::Bool = true

    num_tracers::Int64 = 1
    
    initial_condition::Any      
    
    # -------------------------------------------------------
    # 7. IO Settings
    # -------------------------------------------------------
    output_path::String                   # Base directory for all outputs
    output_filename::String               # Full path for the main NetCDF
    logger::String                        # Full path for the log file
    
    do_raw_output::Bool = true
    pressure_levels::Vector{Float64} = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0]
    vars_to_output::Vector{Symbol} 
    output_interval::Int64

    physics_params::Dict{String, Any} 
end

end