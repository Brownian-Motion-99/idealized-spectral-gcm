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
    base_dict = Dict{Symbol, Tuple{String, String, String, Int}}(
        :vor    => VarMeta("vor",    "s-1",        "Relative Vorticity",         "atmosphere_relative_vorticity",       3),
        :div    => VarMeta("div",    "s-1",        "Divergence",                 "divergence_of_wind",                  3),
        :t      => VarMeta("ta",     "K",          "Air Temperature",            "air_temperature",                     3),
        :ps     => VarMeta("ps",     "Pa",         "Surface Pressure",           "surface_air_pressure",                2),
        :q      => VarMeta("hus",    "1",          "Specific Humidity",          "specific_humidity",                   3),
        :u      => VarMeta("ua",     "m s-1",      "Eastward Wind",              "eastward_wind",                       3),
        :v      => VarMeta("va",     "m s-1",      "Northward Wind",             "northward_wind",                      3),
        :w      => VarMeta("wap",    "Pa s-1",     "Vertical Pressure Velocity", "lagrangian_tendency_of_air_pressure", 3),
        :z      => VarMeta("zg",     "m2 s-2",     "Geopotential",               "geopotential",                        3),
        :precip => VarMeta("pr",     "kg m-2 s-1", "Precipitation Rate",         "precipitation_flux",                  2),
        :shflx  => VarMeta("hfss",   "W m-2",      "Surface Sensible Heat Flux", "surface_upward_sensible_heat_flux",   2),
        :lhflx  => VarMeta("hfls",   "W m-2",      "Surface Latent Heat Flux",   "surface_upward_latent_heat_flux",     2),
        # :vor    => ("vor",     "1/s",     "Relative Vorticity",             3),
        # :div    => ("div",     "1/s",     "Divergence",                     3),
        # :t      => ("t",       "K",       "Temperature",                    3),
        # :ps     => ("ps",      "Pa",      "Surface Pressure",               2),
        # :q      => ("q",       "kg/kg",   "Specific Humidity",              3),
        # :u      => ("u",       "m/s",     "Zonal Wind",                     3),
        # :v      => ("v",       "m/s",     "Meridional Wind",                3),
        # :w      => ("w",       "Pa/s",    "Vertical Pressure Velocity",     3),
        # :p      => ("p",       "Pa",      "Pressure",                       3),
        # :z      => ("z",       "m^2/s^2", "Geopotential",                   3),
        # :lnps   => ("lnps",    "numeric", "Log Surface Pressure",           2),
        # :t_eq   => ("t_eq",    "K",       "Equilibrium Temperature",        3),
        # :shflx  => ("shflx",   "W/m^2",   "Sensible Heat Flux",             2),
        # :lhflx  => ("lhflx",   "W/m^2",   "Latent Heat Flux",               2),
        # :precip => ("precip",  "mm",      "Pseudo-adiabatic Precipitation", 2),
        # Tendencies...
        :du     => ("du_dt",   "m/s^2",   "Zonal Wind Tendency",            3),
        :dv     => ("dv_dt",   "m/s^2",   "Meridional Wind Tendency",       3),
        :dt     => ("dt_dt",   "K/s",     "Temperature Tendency",           3),
        :dvor   => ("dvor_dt", "1/s^2",   "Relative Vorticity Tendency",    3),
        :ddiv   => ("ddiv_dt", "1/s^2",   "Divergence Tendency",            3),
        :dps    => ("dps_dt",  "K/s",     "Surface Pressure Tendency",      2),
        :dq     => ("dq_dt",   "kg/kg/s", "Specific Humidity Tendency",     3)
    )

    # Additional passive tracers
    for i in 1:10
        sym = Symbol("tr$i")
        base_dict[sym] = ("tr$i", "kg/kg", "Passive Tracer $i", 3)
    end

    return base_dict
end

Get_Var_Info() = Get_Var_Info(Val(:PrimitiveEquation))



# --- Barotropic Mode ---
function Get_Var_Info(::Val{:Barotropic})
    return Dict{Symbol, Tuple{String, String, String, Int}}(
        :vor  => ("vor",  "1/s",     "Relative Vorticity", 2),
        :u    => ("u",    "m/s",     "Zonal Wind",         2),
        :v    => ("v",    "m/s",     "Meridional Wind",    2),
        :ke   => ("ke",   "m^2/s^2", "Kinetic Energy",     2),
        :dvor => ("dvor", "1/s^2",   "Vorticity Tendency", 2)
    )
end



# --- Shallow Water Mode ---
function Get_Var_Info(::Val{:ShallowWater})
    return Dict{Symbol, Tuple{String, String, String, Int}}(
        :h    => ("h",    "m",       "Geopotential Height", 2),
        :u    => ("u",    "m/s",     "Zonal Wind",          2),
        :v    => ("v",    "m/s",     "Meridional Wind",     2),
        :vor  => ("vor",  "1/s",     "Relative Vorticity",  2),
        :div  => ("div",  "1/s",     "Divergence",          2),
        :pv   => ("pv",   "s/m",     "Potential Vorticity", 2),
        :dh   => ("dh_dt","m/s",     "Height Tendency",     2)
    )
end

end