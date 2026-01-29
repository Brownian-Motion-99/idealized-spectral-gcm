module Vert_Coordinate_Module

import Base
using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module

export Vert_Coordinate, Compute_Vert_Coord, Vert_Advection!, Mass_Weighted_Global_Integral, qv_Global_Integral_to_e

"""
There are nd levels (nd+1 interfaces)

a and b should be dimensioned by the number of interfaces = 1 + nd

At these interfaces 
pk = ak + bk*p_surf

here ak is actually ak*p_ref
where p_ref  is a constant reference pressure and p_surf is the instantaneous surface pressure
"""

mutable struct Vert_Coordinate
    nλ::Int64
    nθ::Int64
    nd::Int64
    vert_coord_option::String
    vert_difference_option::String
    vert_advect_scheme::String 
    
    
    p_ref::Float64
    zero_top::Bool
    
    # ak[nd+1] = 0, bk[nd+1] = 1, bk[1] = 0
    ak::Array{Float64, 1}
    bk::Array{Float64, 1}
    
    Δak::Array{Float64, 1}
    Δbk::Array{Float64, 1}
    
    # memory container
    flux::Array{ComplexF64, 3}
    vert_integral::Array{Float64, 3}
        
end



function Vert_Coordinate(
    nλ::Int64, nθ::Int64, nd::Int64,
    vert_coord_option::String, vert_difference_option::String, vert_advect_scheme::String,
    p_ref::Float64 = 101325., zero_top::Bool = true,
    scale_heights::Float64 = 4.0, surf_res::Float64 = 1.0, 
    p_press::Float64 = 0.1,  p_sigma::Float64 = 0.3,  exponent::Float64 = 2.5
)
    
    ak, bk = Compute_Vert_Coord(nd, vert_coord_option, p_ref, zero_top, scale_heights, surf_res, p_press,  p_sigma, exponent)
    Δak, Δbk = ak[2:nd+1]-ak[1:nd], bk[2:nd+1]-bk[1:nd]
    
    flux = zeros(Float64, nλ, nθ, nd+1)
    vert_integral = zeros(Float64, nλ, nθ, 1)
    Vert_Coordinate(nλ, nθ, nd, vert_coord_option, vert_difference_option, vert_advect_scheme, p_ref, zero_top, ak, bk, Δak, Δbk, 
    flux, vert_integral)
    
end



function Compute_Vert_Coord(
    nd::Int64, vert_coord_option::String,
    p_ref::Float64 = 101510., 
    zero_top::Bool = true,
    scale_heights::Float64 = 4.0, surf_res::Float64 = 1.0, 
    p_press::Float64 = 0.1,  p_sigma::Float64 = 0.3,  exponent::Float64 = 2.5
)
    
    
    if (vert_coord_option == "even_sigma") 
        a, b = Compute_Even_Sigma(nd)

    elseif (vert_coord_option == "uneven_sigma") 
        a, b = Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, true)
    
    elseif (vert_coord_option == "simmons_and_burridge")
        a, b = Base.invokelatest(Compute_Simmons_Burridge, nd)
    
    elseif (vert_coord_option == "hybrid") 
        a_sigma, b_sigma = Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, false)
        b_press, a_press = Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, false)
        trans = Transition(nd, b_sigma, p_sigma, p_press)
        a = p_ref * (a_sigma.*trans + a_press.*(1.0 .- trans))
        b = b_sigma.*trans + b_press.*(1.0 .- trans)
    
    elseif (vert_coord_option == "mcm") 
        a, b = Compute_Old_Model_Sigma()
    
    elseif (vert_coord_option == "v197") 
        a, b = Compute_V197_Sigma()
        
    else
        error("Compute_Vert_Coord ", vert_coord_option, "is not a valid value for option")
    end 
    
    return a, b
end 



function Transition(nd::Float64, p::Array{Float64, 1}, p_sigma::Float64, p_press::Float64) 
    
    trans = zeros(Float64, nd+1)
    
    for k = 1:nd+1
        if (p[k] <= p_press)
            trans[k] = 0.0
        elseif (p[k] >= p_sigma) 
            trans[k] = 1.0
        else
            x  = p[k]    - p_press
            xx = p_sigma - p_press
            trans[k] = (sin(0.5*pi*x/xx))^2
        end
    end
    
    return trans
end 



function Compute_Even_Sigma(nd::Int64)
    
    a = zeros(Float64, nd+1)
    b = Array(LinRange(0, 1.0, nd+1))
    return a, b
end 



function Compute_Uneven_Sigma(nd::Int64,  scale_heights::Float64, surf_res::Float64, exponent::Float64, zero_top::Bool)
    """
    ζ = 1 - (k-1]/nd
    b = exp( -(surf_res*ζ + (1.0 - surf_res)*(ζ^exponent))* scale_heights)
    """
    
    a = zeros(Float64, nd+1)
    
    # ζ goes from 1.0 (Top) to 0.0 (Surface)
    ζ = Array(LinRange(1.0, 0.0, nd+1))
    
    # z goes from -1.0 (Top) to 0.0 (Surface)
    z = -(surf_res*ζ + (1.0 - surf_res)*(ζ.^exponent))
    
    # --- ERROR WAS HERE ---
    # Old: b = exp.(-z*scale_heights)  -> exp(positive) -> >1.0
    # New: b = exp.(z*scale_heights)   -> exp(negative) -> <1.0
    b = exp.(z*scale_heights)
    # ----------------------
    
    b[nd+1] = 1.0
    
    if (zero_top) 
        b[1] = 0.0
    end
    return a, b
    
end 



function Compute_V197_Sigma()
    
    nd = 18 
    
    a = zeros(Float64, nd+1)
    b = [0.0; 0.0089163; 0.0342936; 0.0740741; 0.1262002; 0.1886145; 0.2592592;  
    0.3360768; 0.4170096; 0.5000000; 0.5829904; 0.6639231; 0.7407407;  
    0.8113854; 0.8737997; 0.9259259; 0.9657064; 0.9910837; 1.0]
    
    
    return a, b
end 



function Compute_Old_Model_Sigma()
    nd = 14
    
    a = zeros(Float64, nd+1)
    b = [0.0; 0.03; 0.0707; 0.1311; 0.2102; 0.3036; 0.4062; 0.5138; 0.6226; 0.7284; 
    0.8255; 0.9066; 0.9640; 0.9933; 1.0]
    
    return a, b
end 



function Compute_Simmons_Burridge(nd::Int64)

    if nd != 20
        @warn "Simmons & Burridge grid is defined for 20 levels. Interpolation or truncation may occur if nd != 20."
    end

    # Coefficients defined at interfaces (nd+1)
    # These are the standard "L20" values used in CCM3/CAM3
    # A_k is in Pascals, B_k is dimensionless.
    
    ak = [
        2.19422, 4.89520, 9.88241, 18.05201, 29.83724,
        44.62333, 61.60586, 78.51243, 77.31270, 75.90131,
        74.24086, 72.28743, 69.98933, 67.28574, 64.10509,
        60.36321, 55.96111, 50.78224, 44.68920, 37.52196,
        0.0
    ] .* 100.0 # Convert hPa to Pa

    bk = [
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.01505, 0.03276,
        0.05359, 0.07810, 0.10695, 0.14088, 0.18080,
        0.22777, 0.28302, 0.34801, 0.42446, 0.51439,
        1.0
    ]
    
    return ak, bk
end



function Vert_Advection!(vert_coord::Vert_Coordinate, r::AbstractArray{Float64,3}, dz::Array{Float64, 3},  w::Array{Float64,3}, Δt::Int64, vert_advect_scheme::String, rdt::Array{Float64,3})
    """
    Consider the coordinate from atmosphere top to the surface,
    top -1---2---3---4-----> bottom
    
    the velocity w is the downward velocity, dz is also from top to the surface
    the advection flux enters each cell (right hand side) is 
    rdt = -w∂r/∂z = -∂wr/∂z + r∂w/∂z 
    = ( [wr]_{k-1/2} -[wr]_{k+1/2})/Δz_k + r_k(w_{k+1/2} - w_{k-1/2})/Δz_k
    
    w is nλ, nθ, nd+1
    r and dz are nλ, nθ, nd
    """
    nd   = vert_coord.nd
    flux = vert_coord.flux
    
    # no flux boundary condition
    flux[:, :, 1]    .= 0.0
    flux[:, :, nd+1] .= 0.0
    
    # 2nd-order centered scheme assuming variable grid spacing ------
    if vert_advect_scheme == "second_centered_wts"
        
        @views begin
            f_in  = flux[:, :, 2:nd]
            w_in  = w[:, :, 2:nd]
            r_up  = r[:, :, 1:nd-1]
            r_dn  = r[:, :, 2:nd]
            dz_up = dz[:, :, 1:nd-1]
            dz_dn = dz[:, :, 2:nd]
            
            @. f_in = w_in * (r_up + (r_dn - r_up)* dz_up / (dz_up + dz_dn))
        end
        
    #  2nd-order centered scheme assuming uniform grid spacing ------
    elseif vert_advect_scheme == "second_centered"
        
        @views begin
            f_in = flux[:, :, 2:nd]
            w_in = w[:, :, 2:nd]
            r_up = r[:, :, 1:nd-1]
            r_dn = r[:, :, 2:nd]

            @. f_in = w_in * (r_up + r_dn) * 0.5
        end
        
    else 
        error("vert_advect_scheme ", vert_advect_scheme, " is not a valid value for option")
    end
    
    @views begin
        rdt_k    = rdt[:, :, 1:nd]
        flux_k   = flux[:, :, 1:nd]
        flux_kp1 = flux[:, :, 2:nd+1]
        r_k      = r[:, :, 1:nd]
        w_k      = w[:, :, 1:nd]
        w_kp1    = w[:, :, 2:nd+1]
        dz_k     = dz[:, :, 1:nd]

        @. rdt_k = (flux_k - flux_kp1 + r_k * (w_kp1 - w_k)) / dz_k
    end
    
end




function Mass_Weighted_Global_Integral(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_data::Array{Float64, 3}, grid_ps::Array{Float64, 3})
    """
    !  This function returns the mass weighted vertical integral of field,
    !  area averaged over the globe. The units of the result are:
    !  (units of field)*(Kg/meters**2)
    """
    
    grav = atmo_data.grav
    nd = vert_coord.nd
    Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
    vert_integral = vert_coord.vert_integral
    
    Δp = similar(grid_ps)
    
    vert_integral .= 0.0
    for k=1:nd
        Δp .= Δak[k] .+ Δbk[k] * grid_ps
        vert_integral .+= grid_data[:,:,k] .* Δp[:,:,1]
    end
    
    mass_weighted_global_integral = Area_Weighted_Global_Mean(mesh, vert_integral)/grav
    
    return mass_weighted_global_integral
end 


function qv_Global_Integral_to_e(vert_coord::Vert_Coordinate, mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data,
    grid_data::Array{Float64, 3}, grid_ps::Array{Float64, 3})
    # The result would be global integral water vapor (e) and the unit is Pa.
    
    grav = atmo_data.grav
    nd = vert_coord.nd
    Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
    vert_integral = vert_coord.vert_integral
    
    Δp = similar(grid_ps)
    
    vert_integral .= 0.0
    for k=1:nd
        Δp .= Δak[k] .+ Δbk[k] * grid_ps
        vert_integral .+= grid_data[:,:,k] .* Δp[:,:,1]
    end
    
    # precipitable water
    Pv = vert_integral ./ grav

    total_water = Area_Weighted_Global_Mean(mesh, Pv) * 4. * pi * atmo_data.radius^2
    return total_water
end 

end