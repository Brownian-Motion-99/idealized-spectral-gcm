using ...Experiment_Configuration
using ...Spectral_Spherical_Mesh_Module
using ...Vert_Coordinate_Module
using ...Atmo_Data_Module
using ...Dyn_Data_Module
using ...Semi_Implicit_Module

function Spectral_Physics!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh,
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data,
    dyn_data::Dyn_Data,
    semi_implicit::Semi_Implicit_Solver,
    physics_params::Dict{String,Any},
)

    grid_δu, grid_δv, grid_δps, grid_δt =
        dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δt
    grid_u_p, grid_v_p, grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_t_p
    grid_p_half, grid_p_full = dyn_data.grid_p_half, dyn_data.grid_p_full
    grid_t_eq = dyn_data.grid_t_eq
    grid_lrf_tendency = dyn_data.grid_lrf_tendency
    grid_bm_t_tendency = dyn_data.grid_bm_t_tendency
    grid_bm_q_tendency = dyn_data.grid_bm_q_tendency
    grid_bm_precip = dyn_data.grid_bm_precip

    grid_δq = dyn_data.grid_δq
    spe_δq = dyn_data.spe_δq

    # Initialize all additive physics tendencies and diagnostics once per call.
    grid_δu .= 0.0
    grid_δv .= 0.0
    grid_δps .= 0.0
    grid_δt .= 0.0
    spe_δq .= 0.0
    grid_δq .= 0.0
    grid_lrf_tendency .= 0.0
    grid_bm_t_tendency .= 0.0
    grid_bm_q_tendency .= 0.0
    grid_bm_precip .= 0.0

    #####################################################################################################
    # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps =
        dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor =
        dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv =
        dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt =
        dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt

    # grid quantities
    grid_u_p, grid_u, grid_u_n = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n =
        dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t_c, grid_t_n = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n

    # related quanties
    grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half,
    dyn_data.grid_lnp_half,
    dyn_data.grid_p_full,
    dyn_data.grid_lnp_full
    grid_dλ_ps, grid_dθ_ps = dyn_data.grid_dλ_ps, dyn_data.grid_dθ_ps
    grid_lnps = dyn_data.grid_lnps

    grid_div, grid_absvor, grid_vor =
        dyn_data.grid_div, dyn_data.grid_absvor, dyn_data.grid_vor
    grid_w_full, grid_M_half = dyn_data.grid_w_full, dyn_data.grid_M_half
    grid_geopots, grid_geopot_full, grid_geopot_half =
        dyn_data.grid_geopots, dyn_data.grid_geopot_full, dyn_data.grid_geopot_half

    grid_energy_full, spe_energy = dyn_data.grid_energy_full, dyn_data.spe_energy

    # moisture pre-process
    spe_q_n = dyn_data.spe_q_n
    spe_q_c = dyn_data.spe_q_c
    spe_q_p = dyn_data.spe_q_p

    grid_q_n = dyn_data.grid_q_n
    grid_q_c = dyn_data.grid_q_c
    grid_q_p = dyn_data.grid_q_p

    spe_δq = dyn_data.spe_δq
    grid_δq = dyn_data.grid_δq

    grid_z_full = dyn_data.grid_z_full
    grid_z_half = dyn_data.grid_z_half

    grid_w_full = dyn_data.grid_w_full

    integrator = semi_implicit.integrator
    Δt = Get_Δt(integrator)

    grid_shflx = dyn_data.grid_shflx
    grid_lhflx = dyn_data.grid_lhflx

    grid_liquid_water_content = dyn_data.grid_liquid_water_content
    grid_precip = dyn_data.grid_precip
    grid_precip .= 0.0
    grid_liquid_water_content .= 0.0

    grid_z_full = dyn_data.grid_z_full
    grid_z_half = dyn_data.grid_z_half
    grid_δq = dyn_data.grid_δq

    K_E = dyn_data.K_E
    ###############################################################################
    # pressure difference
    grid_Δp = dyn_data.grid_Δp
    # temporary variables
    grid_δQ = dyn_data.grid_d_full1

    # Surface wind speed, surface geopotential, density
    V_c, za, rho = Calculate_V_c_za_rho(
        atmo_data,
        grid_p_half,
        grid_p_full,
        grid_ps,
        grid_u,
        grid_v,
        grid_t_c,
        grid_q_c,
    )

    # Surface sensible heat fluxes
    if physics_params["do_Sensible_Heating"]
        C_H = physics_params["C_H"]::Float64
        Sensible_Heating!(mesh, atmo_data, grid_t_c, grid_shflx, V_c, za, rho, Δt, C_H)
        Trans_Grid_To_Spherical!(mesh, grid_t_c, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t_c)
    end

    # Surface latent heat fluxes
    if physics_params["do_Surface_Evaporation"] && config.moisture_processes
        C_E = physics_params["C_E"]::Float64
        Surface_Evaporation!(
            mesh,
            atmo_data,
            grid_ps,
            grid_q_c,
            grid_lhflx,
            V_c,
            za,
            rho,
            Δt,
            C_E,
        )
        Trans_Grid_To_Spherical!(mesh, grid_q_c, spe_q_c)
        Trans_Spherical_To_Grid!(mesh, spe_q_c, grid_q_c)
    end

    # PBL mixing for temperature and moisture
    if physics_params["do_Implicit_PBL_Scheme"]
        C_D = physics_params["C_D"]::Float64
        Implicit_PBL_Mixing!(
            atmo_data,
            grid_p_full,
            grid_p_half,
            grid_t_c,
            grid_q_c,
            K_E,
            V_c,
            za,
            rho,
            physics_params,
            Δt,
            C_D,
        )

        Trans_Grid_To_Spherical!(mesh, grid_t_c, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t_c)

        Trans_Grid_To_Spherical!(mesh, grid_q_c, spe_q_c)
        Trans_Spherical_To_Grid!(mesh, spe_q_c, grid_q_c)
    end

    # Held-Suarez
    if get(physics_params, "do_HS_Forcing", false)
        HS_Forcing!(
            atmo_data,
            Δt,
            86400,
            mesh.sinθ,
            grid_u_p,
            grid_v_p,
            grid_p_half,
            grid_p_full,
            grid_t_p,
            grid_δu,
            grid_δv,
            grid_t_eq,
            grid_δt,
            physics_params,
        )
    end

    # Betts-Miller first, followed by large-scale condensation 
    # diagnosed from the temporary post-convection state.
    do_betts_miller = get(physics_params, "do_Betts_Miller", false)
    do_lscale_cond = get(physics_params, "do_Lscale_Cond", false)
    if do_betts_miller || do_lscale_cond
        config.moisture_processes ||
            error("moist convection and condensation require moisture_processes = true")

        if do_betts_miller
            haskey(physics_params, "BM_state") ||
                error("BM_state was not initialized before time integration")
            bm_state = physics_params["BM_state"]::Betts_Miller_State
            Δt <= bm_state.tau || throw(
                ArgumentError(
                    "Betts-Miller requires effective_dt <= bm_tau; " *
                    "got $Δt s and $(bm_state.tau) s",
                ),
            )
        else
            bm_state = nothing
        end
        if do_lscale_cond
            haskey(physics_params, "L") ||
                error("large-scale condensation requires physics_params[\"L\"]")
            heating_scale = physics_params["L"]
        else
            heating_scale = 1.0
        end

        Moist_Physics!(
            do_betts_miller,
            do_lscale_cond,
            bm_state,
            heating_scale,
            atmo_data,
            grid_t_c,
            grid_q_c,
            grid_p_full,
            grid_p_half,
            Δt,
            grid_bm_t_tendency,
            grid_bm_q_tendency,
            grid_bm_precip,
            dyn_data.grid_d_full1,
            dyn_data.grid_d_full2,
            grid_liquid_water_content,
            grid_δt,
            grid_δq,
            grid_precip,
        )
    end

    # Linear response function for moisture-radiative feedback. 
    # Diagnose it from the same current humidity state used by Betts-Miller.
    if get(physics_params, "do_LRF", false)
        config.moisture_processes || error("LRF requires moisture_processes = true")
        haskey(physics_params, "LRF_state") ||
            error("LRF_state was not initialized before time integration")

        LRF!(
            physics_params["LRF_state"]::LRF_State,
            grid_q_c,
            grid_lrf_tendency,
            config.day_to_sec,
        )
        grid_δt .+= grid_lrf_tendency
    end

end
