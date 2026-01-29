module Spectral_Dynamics_Module

using Base.Threads
using LinearAlgebra
using Printf
using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module
using ..Vert_Coordinate_Module
using ..Dyn_Data_Module
using ..Semi_Implicit_Module
using ..Time_Integrator_Module
using ..Press_And_Geopot_Module
using ..Atmos_Param_Module
using ..Experiment_Configuration

using Statistics
using Interpolations
export Compute_Corrections_Init, Compute_Corrections!, Four_In_One!, Spectral_Dynamics!, Get_Topography!, Spectral_Initialize_Fields!, Spectral_Dynamics_Physics!, Atmosphere_Update!

function Compute_Corrections_Init(
    mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, 
    atmo_data::Atmo_Data,
    grid_ps_p::Array{Float64, 3},
    grid_energy_temp::Array{Float64, 3}, grid_u_p::Array{Float64, 3}, grid_v_p::Array{Float64, 3}, grid_t_p::Array{Float64, 3}, grid_δu::Array{Float64, 3}, grid_δv::Array{Float64, 3}, grid_δt::Array{Float64, 3},  
    grid_q_p::Array{Float64, 3}, grid_δq::Array{Float64,3},
    Δt::Int64,
)
    
    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    cp_air = atmo_data.cp_air

    if (do_mass_correction) 
        mean_ps_p = Area_Weighted_Global_Mean(mesh, grid_ps_p)
    end
    
    if (do_energy_correction) 
        grid_energy_temp .= 0.5 * ((grid_u_p + Float64(Δt) * grid_δu).^2 + (grid_v_p + Float64(Δt) * grid_δv).^2) + cp_air * (grid_t_p + Float64(Δt) * grid_δt)
        mean_energy_p     = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_p)
    end

    if (do_water_correction)
        mean_moisture_p =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_q_p .+ grid_δq * Float64(Δt), grid_ps_p)
    end
    
    return mean_ps_p, mean_energy_p, mean_moisture_p 
end 



function Compute_Corrections!(
    mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, 
    atmo_data::Atmo_Data,
    mean_ps_p::Float64, grid_ps_n::Array{Float64, 3},spe_lnps_n::Array{ComplexF64, 3}, 
    mean_energy_p::Float64, grid_energy_temp::Array{Float64, 3}, grid_u_n::Array{Float64, 3}, grid_v_n::Array{Float64, 3}, grid_t_n::Array{Float64, 3}, spe_t_n::Array{ComplexF64, 3},
    mean_moisture_p::Float64, grid_q_n::Array{Float64, 3}
)

    do_mass_correction, do_energy_correction, do_water_correction = atmo_data.do_mass_correction, atmo_data.do_energy_correction, atmo_data.do_water_correction
    cp_air, grav = atmo_data.cp_air, atmo_data.grav
    
    if (do_mass_correction) 
        mean_ps_n              = Area_Weighted_Global_Mean(mesh, grid_ps_n)
        mass_correction_factor = mean_ps_p / mean_ps_n
        grid_ps_n            .*= mass_correction_factor
        spe_lnps_n[1,1,1]     += log(mass_correction_factor)
    end
    
    if (do_energy_correction) 
        grid_energy_temp      .= 0.5 * (grid_u_n.^2 + grid_v_n.^2) + cp_air * grid_t_n
        mean_energy_n          = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_energy_temp, grid_ps_n)
        
        temperature_correction = grav*(mean_energy_p - mean_energy_n) / (cp_air * mean_ps_p)
        grid_t_n             .+= temperature_correction
        spe_t_n[1,1,:]       .+= temperature_correction
    end

    if (do_water_correction) 
        grid_q_n[grid_q_n .< 0.] .=  0.
        mean_moisture_n           =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_q_n, grid_ps_n)
        grid_q_n                .*=  mean_moisture_p ./ mean_moisture_n 
        mean_moisture_n           =  Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, grid_q_n, grid_ps_n) 
    end
    
end 



"""
compute vertical mass flux and velocity 
grid_M_half[:,:,k+1] = downward mass flux/per unit area across the K+1/2
grid_w_full[:,:,k]   = dp/dt vertical velocity 

update residuals
grid_δps[:,:,k]  += -∑_{r=1}^nd Dr = -∑_{r=1}^nd ∇(vrΔp_r)
grid_δt[:,:,k]   += κTw/p 
(grid_δu[:,:,k], grid_δv[:,:,k]) -= RT ∇p/p 

!  cell boundary. This is the "vertical velocity" in the hybrid coordinate system.
!  When vertical coordinate is pure sigma: grid_M_half = grid_ps*d(sigma)/dt
"""
function Four_In_One!(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
    grid_div::Array{Float64,3}, grid_u::Array{Float64,3}, grid_v::Array{Float64,3}, 
    grid_ps::Array{Float64,3},  grid_Δp::Array{Float64,3}, grid_lnp_half::Array{Float64,3}, grid_lnp_full::Array{Float64,3}, grid_p_full::Array{Float64,3},
    grid_dλ_ps::Array{Float64,3}, grid_dθ_ps::Array{Float64,3}, 
    grid_t::Array{Float64,3}, 
    grid_M_half::Array{Float64,3}, grid_w_full::Array{Float64,3}, 
    grid_δu::Array{Float64,3}, grid_δv::Array{Float64,3}, grid_δps::Array{Float64,3}, grid_δt::Array{Float64,3}, grid_δq::Array{Float64,3}
)
    
    # Unpack parameters
    rdgas, cp_air          = atmo_data.rdgas, atmo_data.cp_air
    nd, bk                 = vert_coord.nd, vert_coord.bk
    Δak, Δbk               = vert_coord.Δak, vert_coord.Δbk
    vert_difference_option = vert_coord.vert_difference_option
    kappa                  = rdgas / cp_air
    nλ, nθ                 = atmo_data.nλ, atmo_data.nθ
    
    @threads for j = 1:nθ
        for i = 1:nλ

            if (vert_difference_option == "simmons_and_burridge")
                
                # dmean_tot = ∇ ∑_{k=1}^{nd} v_k * Δp_k = ∑_{k=1}^{nd} D_k
                dmean_tot = 0.0

                for k = 1:nd

                    Δp     = grid_Δp[i, j, k]
                    Δlnp_p = grid_lnp_half[i, j, k+1] - grid_lnp_full[i, j, k]
                    Δlnp_m = grid_lnp_full[i, j, k]   - grid_lnp_half[i, j, k]
                    Δlnp   = grid_lnp_half[i, j, k+1] - grid_lnp_half[i, j, k]

                    # Solenoidal
                    # ∇p_k/p = [(lnp_k - lnp_{k-1/2}) * ∇p_{k-1/2} + (lnp_{k+1/2} - lnp_k) * ∇p_{k+1/2}] / Δpk
                    #        = [(lnp_k - lnp_{k-1/2}) * b_{k-1/2} + (lnp_{k+1/2} - lnp_k) * b_{k+1/2}] / Δpk * ∇ps
                    #        = x1 * ∇ps
                    x1      = (bk[k] * Δlnp_m + bk[k+1] * Δlnp_p) / Δp
                    dlnp_dλ = x1 * grid_dλ_ps[i, j, 1]
                    dlnp_dθ = x1 * grid_dθ_ps[i, j, 1]

                    # (grid_δu, grid_δv) -= RT/p * ∇p
                    t_k               = grid_t[i, j, k]
                    grid_δu[i, j, k] -= rdgas * t_k * dlnp_dλ
                    grid_δv[i, j, k] -= rdgas * t_k * dlnp_dθ

                    # dmean = ∇⋅(v_k * Δp_k) = div_k * Δp_k + v_k * (Δbk_k * ∇p_s)
                    dmean = grid_div[i, j, k] * Δp + Δbk[k] * (grid_u[i, j, k] * grid_dλ_ps[i, j, 1] + grid_v[i, j, k] * grid_dθ_ps[i, j, 1])

                    # energy conservation for temperature
                    # w/p = dlnp/dt = ∂lnp/∂t + ∂lnp/∂σ dσ + v⋅∇lnp
                    # ∂ξ_k/∂σ dσ = [M_{k+1/2}(ξ_k+1/2 - ξ_k) + M_{k-1/2}(ξ_k - ξ_k-1/2)]/Δp_k
                    # weight the same way
                    # vertical advection operator (M is the downward speed)
                    # ∂lnp_k/∂σ dσ = [M_{k+1/2} (lnp_k+1/2 - lnp_k) + M_{k-1/2} (lnp_k - lnp_k-1/2)] / Δp_k
                    # ∂lnp/∂t = 1/p ∂p/∂t = [∂p/∂t_{k+1/2} (lnp_k+1/2 - lnp_k) + ∂p/∂t_{k-1/2} (lnp_k - lnp_k-1/2)] / Δp_k
                    # As we know
                    # ∂p/∂t_{k+1/2} = -∑_{r=1}^k Dr - M_{k+1/2}
                    
                    # ∂lnp/∂t + dσ ∂lnp/∂σ =  [(-∑_{r=1}^k Dr)(lnp_k+1/2 - lnp_k) + (-∑_{r=1}^{k-1} Dr)(lnp_k - lnp_k-1/2)]/Δp_k
                    #                      = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k
                    
                    x5 = -(dmean_tot * Δlnp + dmean * Δlnp_p) / Δp + grid_u[i, j, k] * dlnp_dλ + grid_v[i, j, k] * dlnp_dθ
                    
                    # grid_δt += κT w/p
                    grid_δt[i, j, k] += kappa * grid_t[i, j, k] * x5
                    
                    # grid_w_full = w
                    grid_w_full[i, j, k] = x5 * grid_p_full[i, j, k]
                    
                    # update dmean_tot to ∑_{r=1}^k ∇⋅(v_r * Δp_r)
                    dmean_tot += dmean
                    
                    # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
                    grid_M_half[i, j, k+1] = -dmean_tot

                end

            else
                error("vert_difference_option ", vert_difference_option, " is not a valid value for option")
            end

            # ∂ps/∂t = -∑_{r=1}^nd ∇(vrΔp_r) = -dmean_tot
            grid_δps[i, j, 1] -= dmean_tot
            
            for k = 1:nd-1
                # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
                grid_M_half[i, j, k+1] += dmean_tot * bk[k+1]
            end
            
            grid_M_half[i, j, 1]    = 0.0
            grid_M_half[i, j, nd+1] = 0.0

        end
    end

end 



"""
The governing equations are
∂div/∂t = ∇ × (A, B) - ∇^2E := f^d                    
∂lnps/∂t= (-∑_k div_k Δp_k + v_k ∇ Δp_k)/ps := f^p    
∂T/∂t = -(u,v)∇T - dσ∂T∂σ + κTw/p + J:= f^t           
Φ = f^Φ                                               

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 

We have 
δdiv = f^d - I^d + I^d
δlnps = f^p - I^p + I^p
δT = f^t - I^t + I^t

"""
function Spectral_Dynamics!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, 
    atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
    semi_implicit::Semi_Implicit_Solver
)
    
    # spectral equation quantities
    spe_lnps_p, spe_lnps_c, spe_lnps_n, spe_δlnps = dyn_data.spe_lnps_p, dyn_data.spe_lnps_c, dyn_data.spe_lnps_n, dyn_data.spe_δlnps
    spe_vor_p, spe_vor_c, spe_vor_n, spe_δvor     = dyn_data.spe_vor_p, dyn_data.spe_vor_c, dyn_data.spe_vor_n, dyn_data.spe_δvor
    spe_div_p, spe_div_c, spe_div_n, spe_δdiv     = dyn_data.spe_div_p, dyn_data.spe_div_c, dyn_data.spe_div_n, dyn_data.spe_δdiv
    spe_t_p, spe_t_c, spe_t_n, spe_δt             = dyn_data.spe_t_p, dyn_data.spe_t_c, dyn_data.spe_t_n, dyn_data.spe_δt
    spe_tracers_p, spe_tracers_c, spe_tracers_n, spe_δtracers = dyn_data.spe_tracers_p, dyn_data.spe_tracers_c, dyn_data.spe_tracers_n, dyn_data.spe_δtracers

    
    # grid quantities
    grid_u_p, grid_u, grid_u_n    = dyn_data.grid_u_p, dyn_data.grid_u_c, dyn_data.grid_u_n
    grid_v_p, grid_v, grid_v_n    = dyn_data.grid_v_p, dyn_data.grid_v_c, dyn_data.grid_v_n
    grid_ps_p, grid_ps, grid_ps_n = dyn_data.grid_ps_p, dyn_data.grid_ps_c, dyn_data.grid_ps_n
    grid_t_p, grid_t, grid_t_n    = dyn_data.grid_t_p, dyn_data.grid_t_c, dyn_data.grid_t_n
    grid_tracers_p, grid_tracers, grid_tracers_n, grid_δtracers = dyn_data.grid_tracers_p, dyn_data.grid_tracers_c, dyn_data.grid_tracers_n, dyn_data.grid_δtracers


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

    grav              = atmo_data.grav
    integrator        = semi_implicit.integrator
    Δt                = Get_Δt(integrator)

    grid_z_full       = dyn_data.grid_z_full
    grid_z_half       = dyn_data.grid_z_half
    grid_δq     = dyn_data.grid_δq 

    K_E               = dyn_data.K_E
    # pressure difference
    grid_Δp             = dyn_data.grid_Δp
    # temporary variables
    grid_δQ             = dyn_data.grid_d_full1
        
    # incremental quantities
    grid_δu, grid_δv, grid_δps, grid_δlnps, grid_δt = dyn_data.grid_δu, dyn_data.grid_δv, dyn_data.grid_δps, dyn_data.grid_δlnps, dyn_data.grid_δt
    integrator          = semi_implicit.integrator
    Δt                  = Get_Δt(integrator)
    ###############################################################################
    mean_ps_p, mean_energy_p, mean_moisture_p = Compute_Corrections_Init(
        mesh, vert_coord,
        atmo_data,
        grid_ps_p,
        grid_energy_full, grid_u_p, grid_v_p, grid_t_p, grid_δu, grid_δv, grid_δt,
        grid_q_p, grid_δq,
        Δt
    )
    
    # compute pressure based on grid_ps -> grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full 
    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)

    # compute ∇ps = ∇lnps * ps
    Compute_Gradients!(mesh, spe_lnps_c,  grid_dλ_ps, grid_dθ_ps)
    grid_dλ_ps .*= grid_ps
    grid_dθ_ps .*= grid_ps

    # compute grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, 
    # except the contributions from geopotential or vertical advection
    Four_In_One!(vert_coord, atmo_data, grid_div, grid_u, grid_v, grid_ps, 
    grid_Δp, grid_lnp_half, grid_lnp_full, grid_p_full,
    grid_dλ_ps, grid_dθ_ps, 
    grid_t, 
    grid_M_half, grid_w_full, grid_δu, grid_δv, grid_δps, grid_δt, grid_δq)

    Compute_Geopotential!(vert_coord, atmo_data, 
    grid_lnp_half, grid_lnp_full,  
    grid_t, 
    grid_geopots, grid_geopot_full, grid_geopot_half, grid_q_c)

    grid_δlnps .= grid_δps ./ grid_ps
    Trans_Grid_To_Spherical!(mesh, grid_δlnps, spe_δlnps)

    # compute vertical advection, todo  finite volume method 
    Vert_Advection!(vert_coord, grid_u, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δu  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_v, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δv  .+= grid_δQ
    Vert_Advection!(vert_coord, grid_t, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
    grid_δt  .+= grid_δQ
    
    # passive tracers
    for i in 1:dyn_data.num_tracers
        @views begin
            spe_p  = spe_tracers_p[:, :, :, i]
            spe_c  = spe_tracers_c[:, :, :, i]
            spe_n  = spe_tracers_n[:, :, :, i]
            spe_δ  = spe_δtracers[:, :, :, i]
            grid_c = grid_tracers[:, :, :, i]
            grid_n = grid_tracers_n[:, :, :, i]
            grid_δ = grid_δtracers[:, :, :, i]
            Vert_Advection!(vert_coord, grid_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme, grid_δQ)
            grid_δ .+= grid_δQ 
            Add_Horizontal_Advection!(mesh, spe_c, grid_u, grid_v, grid_δ) 
            Trans_Grid_To_Spherical!(mesh, grid_δ, spe_δ)
            Compute_Spectral_Damping!(integrator, spe_c, spe_p, spe_δ)
            Filtered_Leapfrog!(integrator, spe_δ, spe_p, spe_c, spe_n)
            Trans_Spherical_To_Grid!(mesh, spe_n, grid_n)
        end
    end
    
    # moisture
    if config.moisture_processes
        Vert_Advection!(vert_coord, grid_q_c, grid_Δp, grid_M_half, Δt, vert_coord.vert_advect_scheme,  grid_δQ)
        grid_δq .+= grid_δQ 
        Add_Horizontal_Advection!(mesh, spe_q_c, grid_u, grid_v, grid_δq) 
        Trans_Grid_To_Spherical!(mesh, grid_δq, spe_δq)
        Compute_Spectral_Damping!(integrator, spe_q_c, spe_q_p, spe_δq)
        Filtered_Leapfrog!(integrator, spe_δq, spe_q_p, spe_q_c, spe_q_n)
        Trans_Spherical_To_Grid!(mesh, spe_q_n, grid_q_n)
    end
    
    # temperature
    Add_Horizontal_Advection!(mesh, spe_t_c, grid_u, grid_v, grid_δt)
    Trans_Grid_To_Spherical!(mesh, grid_δt, spe_δt)
    

    grid_absvor = dyn_data.grid_absvor
    Compute_Abs_Vor!(grid_vor, atmo_data.coriolis, grid_absvor)
    
    grid_δu .+=  grid_absvor .* grid_v
    grid_δv .-=  grid_absvor .* grid_u
    
    Vor_Div_From_Grid_UV!(mesh, grid_δu, grid_δv, spe_δvor, spe_δdiv)

    grid_energy_full .= grid_geopot_full .+ 0.5 * (grid_u.^2 + grid_v.^2)
    Trans_Grid_To_Spherical!(mesh, grid_energy_full, spe_energy)
    Apply_Laplacian!(mesh, spe_energy)
    spe_δdiv .-= spe_energy
    
    
    
    Implicit_Correction!(semi_implicit, vert_coord, atmo_data,
    spe_div_c, spe_div_p, spe_lnps_c, spe_lnps_p, spe_t_c, spe_t_p, 
    spe_δdiv, spe_δlnps, spe_δt)

    
    Compute_Spectral_Damping!(integrator, spe_vor_c, spe_vor_p, spe_δvor)
    Compute_Spectral_Damping!(integrator, spe_div_c, spe_div_p, spe_δdiv)
    Compute_Spectral_Damping!(integrator, spe_t_c, spe_t_p, spe_δt)


    Filtered_Leapfrog!(integrator, spe_δvor, spe_vor_p, spe_vor_c, spe_vor_n)
    Filtered_Leapfrog!(integrator, spe_δdiv, spe_div_p, spe_div_c, spe_div_n)
    Filtered_Leapfrog!(integrator, spe_δlnps, spe_lnps_p, spe_lnps_c, spe_lnps_n)
    Filtered_Leapfrog!(integrator, spe_δt, spe_t_p, spe_t_c, spe_t_n)
    
    
    Trans_Spherical_To_Grid!(mesh, spe_vor_n, grid_vor)
    Trans_Spherical_To_Grid!(mesh, spe_div_n, grid_div)
    UV_Grid_From_Vor_Div!(mesh, spe_vor_n, spe_div_n, grid_u_n, grid_v_n)
    Trans_Spherical_To_Grid!(mesh, spe_lnps_n, grid_lnps)
    grid_ps_n .= exp.(grid_lnps)
    Trans_Spherical_To_Grid!(mesh, spe_t_n, grid_t_n) 


    Compute_Corrections!(
        mesh, vert_coord,
        atmo_data,
        mean_ps_p, grid_ps_n, spe_lnps_n,
        mean_energy_p, grid_energy_full, grid_u_n, grid_v_n, grid_t_n, spe_t_n,
        mean_moisture_p, grid_q_n
    )

    Time_Advance!(dyn_data)

    Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full)
    
    return 
end 



function Get_Topography!(grid_geopots::Array{Float64, 3}, warm_start_file_name::String = "None")
    if warm_start_file_name == "None" # load warm start file
        grid_geopots .= 0.0
    end

    if warm_start_file_name != "None" # load warm start file
        read_file     = load(warm_start_file_name)
        grid_geopots .= read_file["grid_geopots_final"][:,:,:,1]
    end
    
    return
end 



function Spectral_Initialize_Fields!(
    mesh::Spectral_Spherical_Mesh, atmo_data::Atmo_Data, vert_coord::Vert_Coordinate, 
    sea_level_ps_ref::Float64, init_t::Float64, grid_geopots::Array{Float64,3},
    dyn_data::Dyn_Data, warm_start_file_name::String = "None"
)

    if warm_start_file_name != "None" # load warm start file
        spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
        spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
        grid_u, grid_v, grid_ps, grid_t           = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
        grid_u_p, grid_v_p, grid_ps_p, grid_t_p   = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
        
        grid_lnps,  grid_vor, grid_div            = dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
        
        grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
        nλ, nθ, nd                                = mesh.nλ, mesh.nθ, mesh.nd
        
        # moisture pre-process
        grid_t_n          = dyn_data.grid_t_n
        spe_q_c     = dyn_data.spe_q_c
        spe_q_p     = dyn_data.spe_q_p 

        grid_q_n    = dyn_data.grid_q_n
        grid_q_c    = dyn_data.grid_q_c
        grid_q_p    = dyn_data.grid_q_p 

        grid_u_n      = dyn_data.grid_u_n
        grid_v_n      = dyn_data.grid_v_n
        # Tendency 
        grid_δu = dyn_data.grid_δu
        grid_δv = dyn_data.grid_δv

        grid_δq = dyn_data.grid_δq
        ####################################################################
        @info warm_start_file_name # to make sure get the correct warmstart_PR.dat
        read_file      = load(warm_start_file_name)        
        grid_u[:,:,:] .= read_file["grid_u_c_final"][:,:,:,1]
        grid_v[:,:,:] .= read_file["grid_v_c_final"][:,:,:,1]  
        grid_t        .= read_file["grid_t_c_final"][:,:,:,1] 
        
        grid_lnps     .= log.(read_file["grid_ps_c_final"][:,:,1,1])
        grid_ps       .= read_file["grid_ps_c_final"][:,:,1,1]

        # grid_ps -> grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full
        Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full)
        ####################################################################        
        num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
        
        # for current variables
        spe_t_c       .= read_file["spe_t_c_final"][:,:,:,1]         
        spe_vor_c[:,:,:] .= read_file["spe_vor_c_final"][:,:,:,1]
        spe_div_c[:,:,:] .= read_file["spe_div_c_final"][:,:,:,1]
        spe_lnps_c    .= (read_file["spe_lnps_c_final"][:,:,1,1])
        grid_vor .= read_file["grid_vor_final"][:,:,:,1]  
        grid_div .= read_file["grid_div_final"][:,:,:,1] 
        ####################################################################        
        # for past variables
        spe_vor_p   .= read_file["spe_vor_p_final"][:,:,:,1]
        spe_div_p   .= read_file["spe_div_p_final"][:,:,:,1]
        spe_lnps_p  .= read_file["spe_lnps_p_final"][:,:,:,1]
        spe_t_p     .= read_file["spe_t_p_final"][:,:,:,1]

        grid_u_p    .= read_file["grid_u_p_final"][:,:,:,1]
        grid_v_p    .= read_file["grid_v_p_final"][:,:,:,1]
        grid_ps_p   .= read_file["grid_ps_p_final"][:,:,:,1]
        grid_t_p    .= read_file["grid_t_p_final"][:,:,:,1]
        ####################################################################        
        # q
        grid_q_n .= read_file["grid_q_n_final"][:,:,:,1] 
        grid_q_c .= read_file["grid_q_c_final"][:,:,:,1]
        grid_q_p .= read_file["grid_q_p_final"][:,:,:,1]
        
        spe_q_c  .= read_file["spe_q_c_final"][:,:,:,1]
        spe_q_p  .= read_file["spe_q_p_final"][:,:,:,1]
        ####################################################################

    end

    if warm_start_file_name == "None" # then use original start
        spe_vor_c, spe_div_c, spe_lnps_c, spe_t_c = dyn_data.spe_vor_c, dyn_data.spe_div_c, dyn_data.spe_lnps_c, dyn_data.spe_t_c
        spe_vor_p, spe_div_p, spe_lnps_p, spe_t_p = dyn_data.spe_vor_p, dyn_data.spe_div_p, dyn_data.spe_lnps_p, dyn_data.spe_t_p
        grid_u, grid_v, grid_ps, grid_t = dyn_data.grid_u_c, dyn_data.grid_v_c, dyn_data.grid_ps_c, dyn_data.grid_t_c
        grid_u_p, grid_v_p, grid_ps_p, grid_t_p = dyn_data.grid_u_p, dyn_data.grid_v_p, dyn_data.grid_ps_p, dyn_data.grid_t_p
        
        grid_lnps,  grid_vor, grid_div =  dyn_data.grid_lnps, dyn_data.grid_vor, dyn_data.grid_div
        
        grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_p_half, dyn_data.grid_Δp, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full
        nλ, nθ, nd = mesh.nλ, mesh.nθ, mesh.nd
                
        # update q
        spe_q_c     = dyn_data.spe_q_c
        spe_q_p     = dyn_data.spe_q_p 
            
        grid_q_c    = dyn_data.grid_q_c
        grid_q_p    = dyn_data.grid_q_p 

        rdgas = atmo_data.rdgas
        grid_t         .=  init_t 

        # dΦ/dlnp = -RT    Δp = -ΔΦ/RT
        grid_lnps[:,:,1] .= log(sea_level_ps_ref) .- grid_geopots[:,:,1] ./ (rdgas * init_t) 
        grid_ps   .= exp.(grid_lnps)
        
        
        spe_div_c .= 0.0
        spe_vor_c .= 0.0
      
        # # initial perturbation
        num_fourier, num_spherical = mesh.num_fourier, mesh.num_spherical
        
        initial_perturbation = 1.0e-7/sqrt(2.0)
        # initial vorticity perturbation used in benchmark code
        # In gfdl spe[i,j] =  myspe[i, i+j-1]*√2
        for k = nd-2:nd
          spe_vor_c[2,5,k] = initial_perturbation
          spe_vor_c[6,9,k] = initial_perturbation
          spe_vor_c[2,4,k] = initial_perturbation  
          spe_vor_c[6,8,k] = initial_perturbation
        end
      
        UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
        
      
        # initial spectral fields (and spectrally-filtered) grid fields
        
        Trans_Grid_To_Spherical!(mesh, grid_t, spe_t_c)
        Trans_Spherical_To_Grid!(mesh, spe_t_c, grid_t)
      
        Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps_c)
        Trans_Spherical_To_Grid!(mesh, spe_lnps_c,  grid_lnps)
        grid_ps .= exp.(grid_lnps)
        
      
        Vor_Div_From_Grid_UV!(mesh, grid_u, grid_v, spe_vor_c, spe_div_c)
      
        UV_Grid_From_Vor_Div!(mesh, spe_vor_c, spe_div_c, grid_u, grid_v)
        
        Trans_Spherical_To_Grid!(mesh, spe_vor_c, grid_vor)
        Trans_Spherical_To_Grid!(mesh, spe_div_c, grid_div)
        
        #update pressure variables for hs forcing
        Pressure_Variables!(vert_coord, grid_ps, grid_p_half, grid_Δp,
        grid_lnp_half, grid_p_full, grid_lnp_full)
        
        
        spe_vor_p  .= spe_vor_c
        spe_div_p  .= spe_div_c
        spe_lnps_p .= spe_lnps_c
        spe_t_p    .= spe_t_c
      
      
        grid_u_p   .= grid_u
        grid_v_p   .= grid_v
        grid_ps_p  .= grid_ps
        grid_t_p   .= grid_t

        # Tracer initialization
        qv0             = 0.018
        θc              = mesh.θc # lat
        phi_hw          = 2 * pi / 9 * deg2rad(40)
        p_hw            = 30000.
        phi             = LinRange(-90, 90, nθ)
        p0              = 100000.
        for k in 1:20
            for j in 1:64
               for i in 1:128
                   grid_q_c[i,j,k] = qv0 * exp(-((grid_p_full[i,j,k]/grid_ps[i,j,1] - 1.)*(p0/p_hw))^2) * exp(-((θc[j])/phi_hw)^4) 
                    
               end            
            end
        end
        grid_q_c[:,:,1] .= 0.
        
        Trans_Grid_To_Spherical!(mesh, grid_q_c, spe_q_c)
        Trans_Spherical_To_Grid!(mesh, spe_q_c, grid_q_c)

        grid_q_p .= grid_q_c
        spe_q_p  .= spe_q_c
    end
     
end 



function Atmosphere_Update!(
    config::Model_Config,
    mesh::Spectral_Spherical_Mesh, vert_coord::Vert_Coordinate, 
    atmo_data::Atmo_Data, dyn_data::Dyn_Data, 
    semi_implicit::Semi_Implicit_Solver,
    physics_params::Dict{String, Any}
)

    Spectral_Physics!(config, mesh, atmo_data, dyn_data, semi_implicit, physics_params)
    Spectral_Dynamics!(config, mesh, vert_coord, atmo_data, dyn_data, semi_implicit)

    grid_ps , grid_Δp, grid_p_half, grid_lnp_half, grid_p_full, grid_lnp_full = dyn_data.grid_ps_c,  dyn_data.grid_Δp, dyn_data.grid_p_half, dyn_data.grid_lnp_half, dyn_data.grid_p_full, dyn_data.grid_lnp_full 
    
    grid_t = dyn_data.grid_t_c
    grid_geopots, grid_z_full, grid_z_half = dyn_data.grid_geopots, dyn_data.grid_z_full, dyn_data.grid_z_half

    grid_q_c = dyn_data.grid_q_c
        
    Compute_Pressures_And_Heights!(atmo_data, vert_coord,     
    grid_ps, grid_geopots, grid_t, 
    grid_p_half, grid_Δp, grid_lnp_half, grid_p_full, grid_lnp_full, grid_z_full, grid_z_half, grid_q_c)

    return
end 

end