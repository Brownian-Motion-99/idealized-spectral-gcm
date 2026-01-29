using ...Experiment_Configuration
using ...Spectral_Spherical_Mesh_Module
using ...Vert_Coordinate_Module
using ...Atmo_Data_Module
using ...Dyn_Data_Module
using ...Semi_Implicit_Module

function Spectral_Physics!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
    semi_implicit::Semi_Implicit_Solver, 
    physics_params::Dict{String, Any}
)
    
    grid_δu, grid_δv, grid_δps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
    grid_u_p, grid_v_p,  grid_t_p       = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_t_p
    grid_p_half, grid_p_full            = dyn_data.grid_p_half, dyn_data.grid_p_full
    grid_t_eq                           = dyn_data.grid_t_eq

    grid_δq                       = dyn_data.grid_δq
    spe_δq                        = dyn_data.spe_δq

    # Initialize ps and q
    grid_δps .= 0.0
    spe_δq   .= 0.
    grid_δq  .= 0.

    #####################################################################################################
     # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor     = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv     = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt             = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    
    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n

    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps                                 = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps                                              = dyn_data.grid_lnps
    
    grid_div, grid_absvor, grid_vor                        = dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half                               = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half       = dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half
    
    grid_energy_full, spe_energy                           = dyn_data.grid_energy_full, dyn_data.spe_energy
    
    # moisture pre-process
    spe_q_n     = dyn_data.spe_q_n
    spe_q_c     = dyn_data.spe_q_c
    spe_q_p     = dyn_data.spe_q_p 
    
    grid_q_n    = dyn_data.grid_q_n
    grid_q_c    = dyn_data.grid_q_c
    grid_q_p    = dyn_data.grid_q_p 
    
    spe_δq      = dyn_data.spe_δq  
    grid_δq     = dyn_data.grid_δq 
    
    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half

    grid_w_full       = dyn_data.grid_w_full

    integrator        = semi_implicit.integrator
    Δt                = Get_Δt(integrator)

    grid_shflx = dyn_data.grid_shflx
    grid_lhflx = dyn_data.grid_lhflx

    grid_liquid_water_content = dyn_data.grid_liquid_water_content
    grid_precip               = dyn_data.grid_precip

    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    grid_δq           = dyn_data.grid_δq 

    K_E               = dyn_data.K_E
    ###############################################################################
    # pressure difference
    grid_Δp = dyn_data.grid_Δp
    # temporary variables
    grid_δQ = dyn_data.grid_d_full1
    
    if config.moisture_processes

        # V_c, za, rho
        V_c, za, rho = Calculate_V_c_za_rho(
            atmo_data,
            grid_p_half, grid_p_full, grid_ps,
            grid_u, grid_v,
            grid_t, grid_q_c
        )
        
        if physics_params["do_Lscale_Cond"]
            grid_precip .= 0.0
            Lscale_Cond!(
                vert_coord,
                atmo_data,
                grid_q_c, grid_δq, grid_liquid_water_content, grid_precip,
                grid_t, grid_δt,
                grid_p_full, grid_ps,
                Δt
            )
            grid_q_c[grid_q_c .< 0] .= 0
        
            grid_q_c .= grid_q_c .- grid_δq .* (2*Δt)
            grid_t   .= grid_t   .+ grid_δt .* (2*Δt)
        
            Trans_Grid_To_Spherical!(mesh, grid_q_c, spe_q_c)
            Trans_Spherical_To_Grid!(mesh, spe_q_c, grid_q_c)
            
            Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
            Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
            
            grid_δq .= 0.
            grid_δt .= 0.
        end

        if physics_params["do_Sensible_Heating"]
            Sensible_Heating!(
                mesh, atmo_data,
                grid_t, grid_shflx,
                V_c, za, rho,
                Δt
            )
            Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
            Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
        end

        if physics_params["do_Surface_Evaporation"]
            Surface_Evaporation!(
                mesh, atmo_data,
                grid_ps,
                grid_q_c, grid_lhflx,
                V_c, za, rho,
                Δt
            )
            Trans_Grid_To_Spherical!(mesh, grid_q_c, spe_q_c)
            Trans_Spherical_To_Grid!(mesh, spe_q_c, grid_q_c)
        end

        if physics_params["do_Implicit_PBL_Scheme"]
            Implicit_PBL_Mixing!(
                atmo_data,
                grid_p_full, grid_p_half,
                grid_t, grid_q_c,
                K_E,
                V_c, za, rho,
                physics_params,
                Δt
            )
        
            Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
            Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
            
            Trans_Grid_To_Spherical!(mesh, grid_q_c, spe_q_c)
            Trans_Spherical_To_Grid!(mesh, spe_q_c, grid_q_c)
        end

    end

    grid_δt .= 0.0
    HS_Forcing!(
        atmo_data, Δt, 86400, mesh.sinθ,
        grid_u_p, grid_v_p,
        grid_p_half, grid_p_full,
        grid_t,
        grid_δu, grid_δv,
        grid_t_eq, grid_δt,
        physics_params
    )

end