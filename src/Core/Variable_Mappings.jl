module Variable_Mappings_Module

using ..Dyn_Data_Module

export Get_Data_Pointer, Get_Dyn_Var_Map

"""
    Get_Data_Pointer(dyn_data, var_name)
    
Returns the raw array pointer. Now supports generic tracers :tr1, :tr2...
"""
function Get_Data_Pointer(dyn_data::Dyn_Data, var_name::Symbol)
    
    # Standard Dynamics
    if var_name == :u; return dyn_data.grid_u_c; end
    if var_name == :v; return dyn_data.grid_v_c; end
    if var_name == :t; return dyn_data.grid_t_c; end
    if var_name == :lnps; return dyn_data.grid_lnps; end
    if var_name == :vor; return dyn_data.grid_vor; end
    if var_name == :div; return dyn_data.grid_div; end
    if var_name == :geopots; return dyn_data.grid_geopots; end
    if var_name == :w; return dyn_data.grid_w_full; end

    # Tendencies
    if var_name == :du; return dyn_data.grid_δu; end
    if var_name == :dv; return dyn_data.grid_δv; end
    if var_name == :dt; return dyn_data.grid_δt; end
    if var_name == :dvor; return dyn_data.grid_δvor; end
    if var_name == :ddiv; return dyn_data.grid_δdiv; end
    if var_name == :dps; return dyn_data.grid_δps; end
    
    # Fluxes
    if var_name == :shflx; return dyn_data.grid_shflx; end
    if var_name == :lhflx; return dyn_data.grid_lhflx; end

    # Tracers (Dynamic Handling)
    # :q is alias for Tracer #1 (Specific Humidity)
    if var_name == :q
        if size(dyn_data.grid_tracers_c, 4) >= 1
            return view(dyn_data.grid_tracers_c, :, :, :, 1)
        else
            return nothing
        end
    end

    if var_name == :dq
        if size(dyn_data.grid_tracers_c, 4) >= 1
            return view(dyn_data.grid_δtracers, :, :, :, 1)
        else
            return nothing
        end
    end

    # Handle :tr1, :tr2, etc.
    s_str = string(var_name)
    if startswith(s_str, "tr")
        try
            idx = parse(Int, s_str[3:end])
            if idx <= size(dyn_data.grid_tracers_c, 4)
                return view(dyn_data.grid_tracers_c, :, :, :, idx)
            end
        catch
            return nothing
        end
    end

    return nothing
end



"""
    Get_Dyn_Var_Map(dyn::Dyn_Data)
    Get_Dyn_Var_Map(dyn::Dyn_Data, ::Val{Model_Type})

Returns a Dictionary mapping symbols (e.g., :u, :h) to the actual data arrays in Dyn_Data.
Used by Initialization (writing) and Output (reading).
"""
# --- Base (Full 3D) ---
function Get_Dyn_Var_Map(dyn::Dyn_Data, ::Val{:PrimitiveEquation})
    
    base_map = Dict{Symbol, Any}(
        :vor    => dyn.grid_vor,
        :div    => dyn.grid_div,
        :t      => dyn.grid_t_c,
        :ps     => dyn.grid_ps_c,
        :u      => dyn.grid_u_c,
        :v      => dyn.grid_v_c,
        :w      => dyn.grid_w_full,
        :p      => dyn.grid_p_full,
        :z      => dyn.grid_geopot_full,
        :lnps   => dyn.grid_lnps,
        :t_eq   => dyn.grid_t_eq,
        :shflx  => dyn.grid_shflx,
        :lhflx  => dyn.grid_lhflx,
        :du     => dyn.grid_δu,
        :dv     => dyn.grid_δv,
        :dt     => dyn.grid_δt,
        :dvor   => dyn.grid_δvor,
        :ddiv   => dyn.grid_δdiv,
        :dps    => dyn.grid_δps,
        :dq     => dyn.grid_δtracers
    )

    # Inject Tracers dynamically
    # Alias :q -> Tracer 1 (if exists)
    if size(dyn.grid_tracers_c, 4) >= 1
        base_map[:q]  = view(dyn.grid_tracers_c, :, :, :, 1)
        base_map[:dq] = view(dyn.grid_δtracers, :, :, :, 1)
    end

    # Inject :tr1, :tr2...
    for i in 1:size(dyn.grid_tracers_c, 4)
        sym = Symbol("tr$i")
        base_map[sym] = view(dyn.grid_tracers_c, :, :, :, i)
    end

    return base_map
end



# --- Barotropic Mode ---
function Get_Dyn_Var_Map(dyn::Dyn_Data, ::Val{:Barotropic})
    # Note: KE calculation creates a temporary array. 
    # Writing to :ke during initialization will NOT update the model state.
    ke = 0.5 .* (dyn.grid_u_c[:,:,1].^2 .+ dyn.grid_v_c[:,:,1].^2)
    
    return Dict{Symbol, Any}(
        :vor  => view(dyn.grid_vor, :, :, 1),
        :u    => view(dyn.grid_u_c, :, :, 1),
        :v    => view(dyn.grid_v_c, :, :, 1),
        :ke   => ke,
        :dvor => view(dyn.grid_δvor, :, :, 1)
    )
end



# --- Shallow Water Mode ---
function Get_Dyn_Var_Map(dyn::Dyn_Data, ::Val{:ShallowWater})
    h_field = view(dyn.grid_ps_c, :, :, 1) # Aliasing grid_ps -> h
    
    # PV is diagnostic (read-only for initialization purposes)
    pv_field = dyn.grid_absvor[:,:,1] ./ h_field

    return Dict{Symbol, Any}(
        :h    => h_field,
        :u    => view(dyn.grid_u_c, :, :, 1),
        :v    => view(dyn.grid_v_c, :, :, 1),
        :vor  => view(dyn.grid_vor, :, :, 1),
        :div  => view(dyn.grid_div, :, :, 1),
        :pv   => pv_field,
        :dh   => view(dyn.grid_δlnps, :, :, 1)
    )
end

end