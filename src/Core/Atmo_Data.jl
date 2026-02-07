module Atmo_Data_Module

export Atmo_Data, Compute_Abs_Vor!

struct Atmo_Data

    name::String

    nλ::Int64
    nθ::Int64
    nd::Int64
    
    do_mass_correction::Bool
    do_energy_correction::Bool
    do_water_correction::Bool
    use_virtual_temperature::Bool

    radius::Float64
    omega::Float64
    grav::Float64
    rdgas::Float64
    kappa::Float64
    rvgas::Float64
    cp_air::Float64
    Lv::Float64
    gamma::Float64
    alpha::Float64

    coriolis::Array{Float64,1}

end



"""
    Atmo_Data(
        name,
        nλ, nθ, nd,
        do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature,
        sinθ,
        radius, omega, grav, rdgas, kappa, rvgas,
        gamma,
        kwargs...
    )

Constructs a new `Atmo_Data` object, initializing configuration parameters.

### Parameters
    - name: Simulation name.
    
    - nλ: Longitudinal grid size.
    - nθ: Latitudinal grid size.
    - nd: Number of vertical layers.
    
    - do_mass_correction: Flag of doing mass correction.
    - do_energy_correction: Flag of doing energy correction.
    - do_water_correction: Flag of doing water correction.
    - use_virtual_temperature: Flag of using virtual temperature.
    
    - sinθ: Sine values of Gaussian grids.

    - radius: Radius of planet.
    - omega: Angular velocity of planet.
    - grav: Gravity.
    - rdgas: Ideal gas constant of dry air.
    - kappa: Poisson constant.
    - rvgas: Ideal gas constant of moist air.
    - gamma: Moist lapse rate.

    - kwargs...

### Returns
    - Atmo_Data: A fully initialized data struct.

### Modified
    - nothing

"""
function Atmo_Data(
    name::String,  
    nλ::Int64, nθ::Int64, nd::Int64, 
    do_mass_correction::Bool, do_energy_correction::Bool, do_water_correction::Bool, use_virtual_temperature::Bool,
    sinθ::Array{Float64,1};
    radius::Float64, omega::Float64 = 7.292e-5, grav::Float64 = 9.80, rdgas::Float64 = 287.04, kappa::Float64 = 2.0/7.0, rvgas::Float64 = 461.50,
    gamma::Float64 = 0.0065,
    kwargs...
)

    coriolis =  2 * omega * sinθ 
    cp_air   = rdgas/kappa
    Lv       = 2.5e6
    alpha    = (rdgas * gamma) / grav

    Atmo_Data(
        name, 
        nλ, nθ, nd,
        do_mass_correction, do_energy_correction, do_water_correction, use_virtual_temperature,
        radius, omega, grav, rdgas, kappa, rvgas, cp_air, Lv, gamma, alpha,
        coriolis
    )
end



"""
    Compute_Abs_Vor!(grid_vor, coriolis, grid_absvor)

Computes absolute vorticity.

### Parameters
    - grid_vor: Gaussian grid vorticity.
    - coriolis: coriolis force at different latitudes.
    - grid_absvor: Gaussian grid absolute vorticity.

### Returns
    - nothing

### Modified
    - grid_absvor

"""
function Compute_Abs_Vor!(grid_vor::Array{Float64, 3}, coriolis::Array{Float64, 1}, grid_absvor::Array{Float64, 3})
    nλ, nθ, nd = size(grid_vor)

    for j = 1:nθ
        grid_absvor[:, j, :] .= grid_vor[:, j, :] .+ coriolis[j]
    end
end

end