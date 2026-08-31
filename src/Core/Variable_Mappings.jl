module Variable_Mappings_Module

using ..Dyn_Data_Module

export Get_Data_Pointer, Get_Dyn_Var_Map

"""
    Get_Data_Pointer(dyn_data, var_name)
    
Returns the raw array pointer.

"""
function Get_Data_Pointer(dyn_data::Dyn_Data, var_name::Symbol)

    # Standard Dynamics
    if var_name == :u
        return dyn_data.grid_u_c
    end
    if var_name == :v
        return dyn_data.grid_v_c
    end
    if var_name == :t
        return dyn_data.grid_t_c
    end
    if var_name == :lnps
        return dyn_data.grid_lnps
    end
    if var_name == :vor
        return dyn_data.grid_vor
    end
    if var_name == :div
        return dyn_data.grid_div
    end
    if var_name == :geopots
        return dyn_data.grid_geopots
    end
    if var_name == :w
        return dyn_data.grid_w_full
    end
    if var_name == :q
        return dyn_data.grid_q_c
    end

    # Tendencies
    if var_name == :du
        return dyn_data.grid_δu
    end
    if var_name == :dv
        return dyn_data.grid_δv
    end
    if var_name == :dt
        return dyn_data.grid_δt
    end
    if var_name == :dvor
        return dyn_data.grid_δvor
    end
    if var_name == :ddiv
        return dyn_data.grid_δdiv
    end
    if var_name == :dps
        return dyn_data.grid_δps
    end
    if var_name == :dq
        return dyn_data.grid_δq
    end
    if var_name == :bm_dt
        return dyn_data.grid_bm_t_tendency
    end
    if var_name == :bm_dq
        return dyn_data.grid_bm_q_tendency
    end

    # Fluxes
    if var_name == :shflx
        return dyn_data.grid_shflx
    end
    if var_name == :lhflx
        return dyn_data.grid_lhflx
    end

    # Pseudo-adiabatic precipitation
    if var_name == :precip
        return dyn_data.grid_precip
    end
    if var_name == :bm_precip
        return dyn_data.grid_bm_precip
    end

    return nothing
end



"""
    Get_Dyn_Var_Map(dyn_data)
    Get_Dyn_Var_Map(dyn_data, ::Val{Model_Type})

Returns a Dictionary mapping symbols (e.g., :u, :h) to the actual data arrays in Dyn_Data.
Used by Initialization (writing) and Output (reading).

"""
# --- Base (Full 3D) ---
function Get_Dyn_Var_Map(dyn_data::Dyn_Data, ::Val{:PrimitiveEquation})

    base_map = Dict{Symbol,Any}(
        :vor => dyn_data.grid_vor,
        :div => dyn_data.grid_div,
        :t => dyn_data.grid_t_c,
        :ps => dyn_data.grid_ps_c,
        :u => dyn_data.grid_u_c,
        :v => dyn_data.grid_v_c,
        :w => dyn_data.grid_w_full,
        :p => dyn_data.grid_p_full,
        :q => dyn_data.grid_q_c,
        :z => dyn_data.grid_geopot_full,
        :lnps => dyn_data.grid_lnps,
        :t_eq => dyn_data.grid_t_eq,
        :lrf_dt => dyn_data.grid_lrf_tendency,
        :bm_dt => dyn_data.grid_bm_t_tendency,
        :bm_dq => dyn_data.grid_bm_q_tendency,
        :bm_precip => dyn_data.grid_bm_precip,
        :shflx => dyn_data.grid_shflx,
        :lhflx => dyn_data.grid_lhflx,
        :precip => dyn_data.grid_precip,
        :du => dyn_data.grid_δu,
        :dv => dyn_data.grid_δv,
        :dt => dyn_data.grid_δt,
        :dvor => dyn_data.grid_δvor,
        :ddiv => dyn_data.grid_δdiv,
        :dps => dyn_data.grid_δps,
        :dq => dyn_data.grid_δq,
    )

    return base_map
end



# --- Barotropic Mode ---
function Get_Dyn_Var_Map(dyn_data::Dyn_Data, ::Val{:Barotropic})
    # Note: KE calculation creates a temporary array. 
    # Writing to :ke during initialization will NOT update the model state.
    ke = 0.5 .* (dyn_data.grid_u_c[:, :, 1] .^ 2 .+ dyn_data.grid_v_c[:, :, 1] .^ 2)

    return Dict{Symbol,Any}(
        :vor => view(dyn_data.grid_vor, :, :, 1),
        :u => view(dyn_data.grid_u_c, :, :, 1),
        :v => view(dyn_data.grid_v_c, :, :, 1),
        :ke => ke,
        :dvor => view(dyn_data.grid_δvor, :, :, 1),
    )
end



# --- Shallow Water Mode ---
function Get_Dyn_Var_Map(dyn_data::Dyn_Data, ::Val{:ShallowWater})
    h_field = view(dyn_data.grid_ps_c, :, :, 1) # Aliasing grid_ps -> h

    # PV is diagnostic (read-only for initialization purposes)
    pv_field = dyn_data.grid_absvor[:, :, 1] ./ h_field

    return Dict{Symbol,Any}(
        :h => h_field,
        :u => view(dyn_data.grid_u_c, :, :, 1),
        :v => view(dyn_data.grid_v_c, :, :, 1),
        :vor => view(dyn_data.grid_vor, :, :, 1),
        :div => view(dyn_data.grid_div, :, :, 1),
        :pv => pv_field,
        :dh => view(dyn_data.grid_δlnps, :, :, 1),
    )
end

end
