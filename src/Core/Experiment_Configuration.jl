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
    model_type::Symbol  # :Barotropic, :Shallow_Water, :Primitive_Equation

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
    # 3. Vertical Coordinate (New Flexibility)
    # -------------------------------------------------------
    # Options: "even_sigma", "uneven_sigma", "hybrid"
    vert_coord_option::Any    
    # Options: "simmons_and_burridge", etc.
    vert_difference_option::Any 
    # Options: "second_centered_wts", etc.
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
    # 5. Composition & Physics
    # -------------------------------------------------------
    L::Float64 = 0.2
    
    num_grid_tracters::Int64 = 1
    num_spe_tracters::Int64  = 1
    
    initial_condition::Any      
    
    output_filename::String
    logger::String
    do_raw_output::Bool = true
    pressure_levels::Vector{Float64} = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0]
    vars_to_output::Vector{Symbol} 
    output_interval::Int64

    physics_params::Dict{String, Any} 
end

end