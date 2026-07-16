module Output_Mappings_Module

export Get_Var_Info

struct VarMeta
    nc_name::String
    units::String
    long_name::String
    std_name::String
    dims::Int
end

"""
    Get_Var_Info()
    Get_Var_Info(::Val{Model_Type})

Returns metadata (Name, Units, Long Name, Dimensions) for variables.
"""
# --- Base (Primitive Equation) ---
function Get_Var_Info(::Val{:PrimitiveEquation})
    base_dict = Dict{Symbol,VarMeta}(
        :vor => VarMeta(
            "vor",
            "s-1",
            "Relative Vorticity",
            "atmosphere_relative_vorticity",
            3,
        ),
        :div => VarMeta("div", "s-1", "Divergence", "divergence_of_wind", 3),
        :t => VarMeta("ta", "K", "Air Temperature", "air_temperature", 3),
        :ps => VarMeta("ps", "Pa", "Surface Pressure", "surface_air_pressure", 2),
        :q => VarMeta("hus", "1", "Specific Humidity", "specific_humidity", 3),
        :u => VarMeta("ua", "m s-1", "Eastward Wind", "eastward_wind", 3),
        :v => VarMeta("va", "m s-1", "Northward Wind", "northward_wind", 3),
        :w => VarMeta(
            "wap",
            "Pa s-1",
            "Vertical Pressure Velocity",
            "lagrangian_tendency_of_air_pressure",
            3,
        ),
        :p => VarMeta("p", "Pa", "Air Pressure", "air_pressure", 3),
        :z => VarMeta("zg", "m2 s-2", "Geopotential", "geopotential", 3),
        :precip =>
            VarMeta("pr", "kg m-2 s-1", "Precipitation Rate", "precipitation_flux", 2),
        :shflx => VarMeta(
            "hfss",
            "W m-2",
            "Surface Sensible Heat Flux",
            "surface_upward_sensible_heat_flux",
            2,
        ),
        :lhflx => VarMeta(
            "hfls",
            "W m-2",
            "Surface Latent Heat Flux",
            "surface_upward_latent_heat_flux",
            2,
        ),
        :t_eq => VarMeta(
            "teq",
            "K",
            "Held-Suarez Equilibrium Temperature",
            "held_suarez_equilibrium_temperature",
            3,
        ),
        :lrf_dt => VarMeta(
            "lrf_dta_dt",
            "K s-1",
            "LRF Temperature Tendency",
            "tendency_of_air_temperature_due_to_linear_response_forcing",
            3,
        ),
        :lnps => VarMeta(
            "lnps",
            "1",
            "Logarithm of Surface Pressure",
            "log_surface_air_pressure",
            2,
        ),
        # Tendencies
        :du => VarMeta(
            "dua_dt",
            "m s-2",
            "Tendency of Eastward Wind",
            "tendency_of_eastward_wind",
            3,
        ),
        :dv => VarMeta(
            "dva_dt",
            "m s-2",
            "Tendency of Northward Wind",
            "tendency_of_northward_wind",
            3,
        ),
        :dt => VarMeta(
            "dta_dt",
            "K s-1",
            "Tendency of Air Temperature",
            "tendency_of_air_temperature",
            3,
        ),
        :dvor => VarMeta(
            "dvor_dt",
            "s-2",
            "Tendency of Relative Vorticity",
            "tendency_of_atmosphere_relative_vorticity",
            3,
        ),
        :ddiv => VarMeta(
            "ddiv_dt",
            "s-2",
            "Tendency of Divergence",
            "tendency_of_divergence_of_wind",
            3,
        ),
        :dps => VarMeta(
            "dps_dt",
            "Pa s-1",
            "Tendency of Surface Pressure",
            "tendency_of_surface_air_pressure",
            2,
        ),
        :dq => VarMeta(
            "dq_dt",
            "s-1",
            "Tendency of Specific Humidity",
            "tendency_of_specific_humidity",
            3,
        ),
    )

    # Additional passive tracers
    for i = 1:10
        sym = Symbol("tr$i")
        base_dict[sym] = VarMeta("tr$i", "1", "Passive Tracer $i", "passive_tracer_$i", 3)
    end

    return base_dict
end

Get_Var_Info() = Get_Var_Info(Val(:PrimitiveEquation))



# --- Barotropic Mode ---
function Get_Var_Info(::Val{:Barotropic})
    return Dict{Symbol,VarMeta}(
        :vor =>
            VarMeta("vor", "s-1", "Relative Vorticity", "atmosphere_relative_vorticity", 2),
        :u => VarMeta("ua", "m s-1", "Eastward Wind", "eastward_wind", 2),
        :v => VarMeta("va", "m s-1", "Northward Wind", "northward_wind", 2),
        :ke => VarMeta("ke", "m2 s-2", "Kinetic Energy", "specific_kinetic_energy", 2),
        :dvor => VarMeta(
            "dvor",
            "s-2",
            "Vorticity Tendency",
            "tendency_of_atmosphere_relative_vorticity",
            2,
        ),
    )
end



# --- Shallow Water Mode ---
function Get_Var_Info(::Val{:ShallowWater})
    return Dict{Symbol,VarMeta}(
        :h => VarMeta("height", "m", "Geopotential Height", "geopotential_height", 2),
        :u => VarMeta("ua", "m s-1", "Eastward Wind", "eastward_wind", 2),
        :v => VarMeta("va", "m s-1", "Northward Wind", "northward_wind", 2),
        :vor =>
            VarMeta("vor", "s-1", "Relative Vorticity", "atmosphere_relative_vorticity", 2),
        :div => VarMeta("div", "s-1", "Divergence", "divergence_of_wind", 2),
        :pv => VarMeta(
            "pv",
            "m-1 s-1",
            "Potential Vorticity",
            "atmosphere_potential_vorticity",
            2,
        ),
        :dh => VarMeta(
            "dh_dt",
            "m s-1",
            "Height Tendency",
            "tendency_of_geopotential_height",
            2,
        ),
    )
end

end
