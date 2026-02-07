module Press_And_Geopot_Module

using Base.Threads
using ..Atmo_Data_Module
using ..Vert_Coordinate_Module

export Compute_Pressures_And_Heights!, Half_Level_Pressures!, Pressure_Variables!, Compute_Geopotential!



"""
    Compute_Pressures_And_Heights!(
        atmo_data, vert_coord,
        grid_ps, grid_geopots, grid_t,
        grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full,
        grid_z_full, grid_z_half, 
        grid_q
    )

Updates the diagnostic pressure and geometric height variables for the entire atmospheric column, 
ensuring hydrostatic balance consistent with the provided vertical coordinate system.

### Parameters
    - atmo_data: Structure containing physical constants for the simulation.
    - vert_coord: Vertical discretization definitions.

    - grid_ps: Surface pressure [nλ, nθ, 1].
    - grid_geopots: Surface geopotential [nλ, nθ, 1].
    - grid_t: Atmospheric temperature [nλ, nθ, nd].

    - grid_p_half: Pressure at interfaces [nλ, nθ, nd+1].
    - grid_Δp: Pressure thickness [nλ, nθ, nd].
    
    - grid_lnp_half: Logarithm pressure at interfaces [nλ, nθ, nd+1].
    - grid_p_full: Pressure at layer centers [nλ, nθ, nd].
    - grid_lnp_full: Logarithm pressure at layer centers [nλ, nθ, nd].

    - grid_z_full: Geopotentail at layer centers [nλ, nθ, nd].
    - grid_z_half: Geopotentail at interfaces [nλ, nθ, nd+1].

    - grid_q: Moisture [nλ, nθ, nd].
    
### Returns
    - nothing
    
### Modified
    - grid_p_half
    - grid_Δp
    - grid_lnp_half
    - grid_p_full
    - grid_lnp_full
    - grid_z_full
    - grid_z_half

"""
function Compute_Pressures_And_Heights!(
    atmo_data::Atmo_Data, vert_coord::Vert_Coordinate,
    grid_ps::Array{Float64, 3}, grid_geopots::Array{Float64, 3}, grid_t::Array{Float64, 3},
    grid_p_half::Array{Float64, 3}, grid_Δp::Array{Float64, 3},
    grid_lnp_half::Array{Float64, 3}, grid_p_full::Array{Float64, 3}, grid_lnp_full::Array{Float64, 3},
    grid_z_full::Array{Float64, 3}, grid_z_half::Array{Float64, 3}, 
    grid_q::Array{Float64, 3}
)

    grav = atmo_data.grav

    Pressure_Variables!(
        vert_coord, 
        grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full
    )

    Compute_Geopotential!(vert_coord, atmo_data, grid_lnp_half, grid_lnp_full, grid_t, grid_geopots, grid_z_full, grid_z_half, grid_q)

    grid_z_full ./= grav
    grid_z_half ./= grav

end



"""
    Half_Level_Pressures!(vert_coord, grid_ps, grid_p_half)

Computes the hydrostatic pressure at vertical layer interfaces (half-levels)

### Parameters
    - vert_coord: Structure defining the vertical discretization.
    - grid_ps: Surface pressure [nλ, nθ, 1].
    - grid_p_half: Pressure at interfaces [nλ, nθ, nd+1].

### Returns
    - nothing
    
### Modified
    - grid_p_half

"""
function Half_Level_Pressures!(vert_coord::Vert_Coordinate, grid_ps::Array{Float64,3}, grid_p_half::Array{Float64,3})

    nd = vert_coord.nd
    bk = vert_coord.bk
    ak = vert_coord.ak

    # pk = ak * pref
    for k = 1:nd+1
        grid_p_half[:, :, k] .= ak[k] .+ bk[k] * grid_ps[:, :, 1]
    end
end



"""
    Pressure_Variables!(
        vert_coord,
        grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full
    )

Computes the full suite of diagnostic pressure variables required for the dynamical core, deriving full-level and interface values from surface pressure.

### Parameters
    - vert_coord: Configuration for the vertical grid.
    
    - grid_ps: Surface pressure [nλ, nθ, 1].
    - grid_p_half: Pressure at interfaces [nλ, nθ, nd+1].
    - grid_Δp: Pressure thickness [nλ, nθ, nd].

    - grid_lnp_half: Logarithm pressure at interfaces [nλ, nθ, nd+1]
    - grid_p_full: Pressure at layer centers [nλ, nθ, nd].
    - grid_lnp_full: Logarithm pressure at layer centers [nλ, nθ, nd]

### Returns
    - nothing

### Modified
    - grid_p_half
    - grid_Δp
    - grid_lnp_half
    - grid_p_full
    - grid_lnp_full

"""
function Pressure_Variables!(
    vert_coord::Vert_Coordinate, 
    grid_ps::Array{Float64,3}, grid_p_half::Array{Float64,3}, grid_Δp::Array{Float64,3},
    grid_lnp_half::Array{Float64,3}, grid_p_full::Array{Float64,3}, grid_lnp_full::Array{Float64,3}
)

    @assert(size(grid_ps)[3] == 1)

    Half_Level_Pressures!(vert_coord, grid_ps, grid_p_half)
    nd = vert_coord.nd
    grid_Δp .= grid_p_half[:, :, 2:nd+1] - grid_p_half[:, :, 1:nd]
    zero_top = vert_coord.zero_top

    if (vert_coord.vert_difference_option == "simmons_and_burridge")

        k_top = (zero_top ? 2 : 1)

        grid_lnp_half[:, :, k_top:nd+1] .= log.(grid_p_half[:, :, k_top:nd+1])

        # lnp_{k} = (p_{k+1/2}lnp_{k+1/2} - p_{k-1/2}lnp_{k-1/2})/Δp_k - 1
        #         = [(p_{k+1/2}-p_{k-1/2})lnp_{k+1/2} + p_{k-1/2}(lnp_{k+1/2} - lnp_{k-1/2})]/Δp_k - 1
        #         = lnp_{k+1/2} + [p_{k-1/2}(lnp_{k+1/2} - lnp_{k-1/2})]/Δp_k - 1
        grid_lnp_full[:, :, k_top:nd] .= grid_lnp_half[:, :, k_top+1:nd+1] .+ grid_p_half[:, :, k_top:nd] .* (grid_lnp_half[:, :, k_top+1:nd+1] - grid_lnp_half[:, :, k_top:nd]) ./ grid_Δp[:, :, k_top:nd] .- 1.0

        if (zero_top)
            grid_lnp_half[:, :, 1] .= 0.0
            grid_lnp_full[:, :, 1] .= grid_lnp_half[:, :, 2] .- 1.0
        end

    else
        error("vert_difference_option ", vert_coord.vert_difference_option, " is not a valid value for option")
    end

    grid_p_full .= exp.(grid_lnp_full)

end



"""
    Compute_Geopotential!(
        vert_coord, atmo_data,
        grid_lnp_half, grid_lnp_full,
        grid_t,
        grid_geopots, grid_geopot_full, grid_geopot_half,
        grid_q
    )

Integrates the hydrostatic relation vertically to compute the geopotential field (Φ) at both layer interfaces and layer centers, 
accounting for moisture effects via virtual temperature.

### Parameters
    - vert_coord: Configuration for the vertical grid.

    - grid_lnp_half: Logarithm pressure at interfaces [nλ, nθ, nd+1]
    - grid_lnp_full: Logarithm pressure at layer centers [nλ, nθ, nd]

    - grid_t: Atmospheric temperature [nλ, nθ, nd].

    - grid_geopots: Surface geopotential [nλ, nθ, 1].
    - grid_geopot_full: Geopotential at layer centers [nλ, nθ, nd].
    - grid_geopots: Geopotential at interfaces [nλ, nθ, nd+1].

    - grid_q: Moisture [nλ, nθ, nd].

### Returns
    - nothing

### Modified
    - grid_geopot_full
    - grid_geopot_half

"""
function Compute_Geopotential!(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    grid_lnp_half::Array{Float64, 3}, grid_lnp_full::Array{Float64, 3},
    grid_t::Array{Float64, 3},
    grid_geopots::Array{Float64, 3}, grid_geopot_full::Array{Float64, 3}, grid_geopot_half::Array{Float64, 3}, 
    grid_q::Array{Float64, 3}
)

    use_virtual_temperature = atmo_data.use_virtual_temperature
    rvgas, rdgas            = atmo_data.rvgas, atmo_data.rdgas
    zero_top                = vert_coord.zero_top
    nd                      = vert_coord.nd

    grid_geopot_half[:, :, nd+1] .= grid_geopots[:, :, 1]

    if zero_top  #todo (pk(1).eq.0.0) then
        k_top = 2
        grid_geopot_half[:, :, 1] .= 0.0
    else
        k_top = 1
    end

    virtual_t = zeros(Float64, size(grid_t))
    if (use_virtual_temperature)
        virtual_t .= grid_t .* (1.0 .+ (rvgas / rdgas - 1.0) .* grid_q)
    else
        virtual_t .= grid_t
    end

    for k = nd:-1:k_top
        #Φ_{k-1/2} = Φ_{k+1/2} + RT_k(ln p_{k+1/2} - ln p_{k-1})
        grid_geopot_half[:, :, k] .= grid_geopot_half[:, :, k+1] .+ rdgas * virtual_t[:, :, k] .* (grid_lnp_half[:, :, k+1] - grid_lnp_half[:, :, k])
    end

    for k = 1:nd
        #Φ_{k} = Φ_{k+1/2} + RT_k(ln p_{k+1/2} - ln p_{k})
        grid_geopot_full[:, :, k] .= grid_geopot_half[:, :, k+1] .+ rdgas * virtual_t[:, :, k] .* (grid_lnp_half[:, :, k+1] - grid_lnp_full[:, :, k])
    end

end

end