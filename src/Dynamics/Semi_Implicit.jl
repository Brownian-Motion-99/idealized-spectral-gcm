module Semi_Implicit_Module

using LinearAlgebra
using ..Vert_Coordinate_Module
using ..Atmo_Data_Module
using ..Time_Integrator_Module
using ..Press_And_Geopot_Module

export Semi_Implicit_Solver, Build_Implicit_Matrices, Build_Wave_Matrices, Update_Init_Step!, Linear_Press_Gradient_δps, Linear_Geopot_δt!,
Linear_Geopot_δps, Linear_Ps_T_δdiv!, Adjust_δlnps_δt_δdiv!, Implicit_Correction!

mutable struct Semi_Implicit_Solver
    
    integrator::Filtered_Leapfrog
    nd::Int64

    # reference pressure 
    ps_ref::Float64
    Δp_ref::Array{Float64,1}
    lnp_half_ref::Array{Float64,1}
    lnp_full_ref::Array{Float64,1}

    # reference temperature
    t_ref::Array{Float64,1}

    # linearized matrices and inverts
    h::Array{Float64,1}
    div_mat::Array{Float64,2}

    num_wavenumbers::Int64
    wave_numbers::Array{Int64,2}
    wave_matrix::Array{Float64,3}

    # spectral memory container
    spe_δdiv_temp::Array{ComplexF64,3}
    spe_δps_temp::Array{ComplexF64,3}
    spe_δt_temp::Array{ComplexF64,3}
    spe_δgeopot_half_temp::Array{ComplexF64,3}
    spe_δgeopot_temp::Array{ComplexF64,3}
    spe_M_half_temp::Array{ComplexF64,3}
    spe_Mdt_half_temp::Array{ComplexF64,3}

end



"""
    Semi_Implicit_Solver(
        vert_coord, atmo_data, integrator, 
        ps_ref, t_ref, wave_numbers
    )

Initializes the semi-implicit solver by linearizing the atmospheric governing equations 
around a reference state and pre-computing the Helmholtz inversion matrices for 
spectral space.

### Parameters
    - vert_coord: Vertical discretization definition.
    - atmo_data: Physical constants and atmospheric parameters.
    - integrator: The Filtered_Leapfrog structure used to determine the time-step scaling (ξ).
    
    - ps_ref: Global constant reference surface pressure (Pa).
    - t_ref: Reference temperature profile [nd] used for linearization.
    - wave_numbers: Matrix of spherical harmonic degrees (n) for the spectral grid.

### Returns
    - Semi_Implicit_Solver: An initialized structure containing pre-computed linearized 
    matrices (h, div_mat), spectral wave matrices, and memory buffers for implicit correction.

"""
function Semi_Implicit_Solver(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, integrator::Filtered_Leapfrog, 
    ps_ref::Float64, t_ref::Array{Float64,1}, wave_numbers::Array{Int64,2}
)

    nd, ak, bk = vert_coord.nd, vert_coord.ak, vert_coord.bk
    # ps_ref_3d  = reshape([ps_ref], (1, 1, 1))
    ps_ref_3d  = fill(ps_ref, (1, 1, 1))
    
    Δp_ref                   = zeros(Float64, (1, 1, nd))
    p_half_ref, lnp_half_ref = zeros(Float64, (1, 1, nd+1)), zeros(Float64, (1, 1, nd+1))
    p_full_ref, lnp_full_ref = zeros(Float64, (1, 1, nd)), zeros(Float64, (1, 1, nd))

    # Reference log-pressures and vertical thicknesses.
    Pressure_Variables!(vert_coord, ps_ref_3d, p_half_ref, Δp_ref, lnp_half_ref, p_full_ref, lnp_full_ref)

    Δp_ref                   = Δp_ref[1, 1, :]
    p_half_ref, lnp_half_ref = p_half_ref[1, 1, :], lnp_half_ref[1, 1, :]
    p_full_ref, lnp_full_ref = p_full_ref[1, 1, :], lnp_full_ref[1, 1, :]

    δlnp_half_ref = zeros(Float64, nd+1)
    δlnp_full_ref = zeros(Float64, nd)
    # derivative d lnp_{k+1/2} dps
    for k = 1:nd+1
        δlnp_half_ref[k] = bk[k] / p_half_ref[k]
    end

    for k = 1:nd
        Δlnp_p = lnp_half_ref[k+1] - lnp_full_ref[k]
        Δlnp_m = lnp_full_ref[k] - lnp_half_ref[k]
        δlnp_full_ref[k] = (bk[k+1] * Δlnp_p + bk[k] * Δlnp_m) / Δp_ref[k]
    end

    # The vertical coupling between divergence, temperature, and surface pressure.
    h, div_mat = Build_Implicit_Matrices(
        vert_coord, atmo_data,
        ps_ref, Δp_ref,
        lnp_half_ref, lnp_full_ref, δlnp_half_ref, δlnp_full_ref, t_ref
    )

    num_fourier, num_spherical = size(wave_numbers) .- 1
    num_wavenumbers            = num_spherical - 1

    # first time step 
    ξ = Get_ξ(integrator)
    laplacian_eigen = integrator.laplacian_eigen
    wave_matrix = Build_Wave_Matrices(num_wavenumbers, laplacian_eigen, div_mat, ξ)


    # spectral memory container
    spe_δdiv_temp         = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, nd)
    spe_δt_temp           = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, nd)
    spe_δps_temp          = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, 1)
    spe_δgeopot_half_temp = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, nd+1)
    spe_δgeopot_temp      = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, nd)
    spe_M_half_temp       = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, nd+1)
    spe_Mdt_half_temp     = Array{ComplexF64,3}(undef, num_fourier+1, num_spherical+1, nd+1)


    Semi_Implicit_Solver(
        integrator, nd, ps_ref, Δp_ref, lnp_half_ref, lnp_full_ref, t_ref,
        h, div_mat, num_wavenumbers, wave_numbers, wave_matrix,
        spe_δdiv_temp, spe_δps_temp, spe_δt_temp, spe_δgeopot_half_temp, spe_δgeopot_temp, spe_M_half_temp, spe_Mdt_half_temp
    )

end



"""
    Build_Implicit_Matrices(
        vert_coord, atmo_data,
        ps_ref, Δp_ref,
        lnp_half_ref, lnp_full_ref, δlnp_half_ref, δlnp_full_ref, t_ref
    )

Linearizes the atmospheric equations to build the vertical operators (τ, γ, ν, h) 
required for the semi-implicit gravity wave correction.

### Parameters
    - vert_coord: Vertical grid configuration.
    - atmo_data: Physical constants (e.g., gas constant R, κ).

    - ps_ref: Reference surface pressure.
    - Δp_ref: Reference layer pressure thicknesses.
    - lnp_half_ref: Reference log-pressure at interfaces.
    - lnp_full_ref: Reference log-pressure at layer centers.
    
    - δlnp_half_ref: Derivative of interface log-pressure with respect to surface pressure.
    - δlnp_full_ref: Derivative of layer center log-pressure with respect to surface pressure.
    
    - t_ref: Reference temperature profile.

### Returns
    - h: Combined pressure gradient and geopotential linearization vector with respect to ln(pₛ).
    - div_mat: The final linearized system matrix (γτ + hνᵀ) representing the vertical 
    structure of the Helmholtz equation.

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 
ν is a row vectors, H1 = RT_ref∂lnp/∂ps, H2 are colume vectors, γ is an upper triangular matrix,  τ is a lower triangular matrix
(1 - ξ^2 ∇^2(γ τ + (H2+H1)ν)δdiv = Δdiv

div_mat = (γ τ + (H2+H1)ν)
wave_matrix[:,:,s] = (I + ξ^2 (γ τ + (H2+H1)ν) s(s+1)/r^2)^{-1}

"""
function Build_Implicit_Matrices(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    ps_ref::Float64, 
    Δp_ref::Array{Float64,1}, lnp_half_ref::Array{Float64,1}, lnp_full_ref::Array{Float64,1},
    δlnp_half_ref::Array{Float64,1}, δlnp_full_ref::Array{Float64,1},
    t_ref::Array{Float64,1}
)

    nd      = vert_coord.nd
    tau     = zeros(Float64, nd, nd)    # temperature with respect to div (adiabatic heating)
    gamma   = zeros(Float64, nd, nd)    # geopotential with respect to T (hydrostatic)
    nu      = zeros(Float64, nd)        # lnps with respect to div (continuity)
    h1      = zeros(Float64, nd)        # pressure gradient with respect to ps
    h2      = zeros(Float64, nd)        # geopotential with respect to ps
    h       = zeros(Float64, nd)        # h1 + h2
    div_mat = zeros(Float64, nd, nd)    # gamma * tau + h nu^T   

    spe_zero_half, spe_zero_full = zeros(ComplexF64, 1, 1, nd+1), zeros(ComplexF64, 1, 1, nd)
    spe_input = zeros(ComplexF64, 1, 1, nd)

    spe_δps, spe_δt, spe_δgeopot_half, spe_δgeopot = zeros(ComplexF64, 1, 1, 1), zeros(ComplexF64, 1, 1, nd), zeros(ComplexF64, 1, 1, nd+1), zeros(ComplexF64, 1, 1, nd)
    δlnp_half_zero, δlnp_full_zero = zeros(Float64, 1, 1, nd + 1), zeros(Float64, 1, 1, nd)

    spe_M_half_temp, spe_Mdt_half_temp = zeros(ComplexF64, 1, 1, nd+1), zeros(ComplexF64, 1, 1, nd+1)
    for k = 1:nd
        spe_input .= 0.0
        spe_input[1, 1, k] = 1.0

        Linear_Ps_T_δdiv!(
            vert_coord, atmo_data, 
            spe_input, 
            ps_ref, 
            Δp_ref, lnp_half_ref, lnp_full_ref, 
            t_ref,
            spe_M_half_temp, spe_Mdt_half_temp, spe_δps, spe_δt
        )

        # -d_lnps/d_div
        nu[k] = -spe_δps[1, 1, 1]

        # -d_T/d_div
        tau[:, k] = -spe_δt[1, 1, :]

        Linear_Geopot_δt!(vert_coord, atmo_data, lnp_half_ref, lnp_full_ref, spe_input, t_ref, spe_δgeopot_half, spe_δgeopot)
        gamma[:, k] = spe_δgeopot[1, 1, :]
    end

    # used only once
    h1 = Linear_Press_Gradient_δps(vert_coord, atmo_data, ps_ref, Δp_ref, lnp_half_ref, lnp_full_ref, t_ref)
    h2 = Linear_Geopot_δps(vert_coord, atmo_data, δlnp_half_ref, δlnp_full_ref, lnp_half_ref, lnp_full_ref, t_ref)

    h = h1 + h2

    div_mat = h * nu' + gamma * tau

    return h, div_mat
end



"""
    Build_Wave_Matrices(num_total_wavenumbers, eigen, div_mat, ξ)

Computes the inverse of the vertical coupling matrix for each spectral total wavenumber L. 
This inversion solves the Helmholtz equation arising from the semi-implicit discretization.

### Parameters
    - num_total_wavenumbers: The maximum total wavenumber (truncation limit).
    - eigen: Eigenvalues (λ) of the Laplacian operator for each wavenumber.
    - div_mat: The vertical structure matrix (B) combining thermal and pressure effects.
    - ξ: The time-integration scaling factor (αΔt).

### Returns
    - wave_matrix: A 3D array [nd, nd, num_total_wavenumbers + 1] containing the 
    inverse matrix for each wavenumber.

"""
function Build_Wave_Matrices(num_total_wavenumbers::Int64, eigen::Array{Float64,2}, div_mat::Array{Float64,2}, ξ::Float64)
    
    nd = size(div_mat)[1]
    wave_matrix = zeros(Float64, nd, nd, num_total_wavenumbers+1)

    for i = 0:num_total_wavenumbers
        factor = ξ^2 * eigen[1, i+1]
        for k = 1:nd
            wave_matrix[k, k, i+1] = 1.0
        end
        wave_matrix[:, :, i+1] .-= factor * div_mat
        wave_matrix[:, :, i+1]  .= inv(wave_matrix[:, :, i+1])
    end

    return wave_matrix
end



function Update_Init_Step!(semi_implicit::Semi_Implicit_Solver)
    integrator = semi_implicit.integrator
    Time_Integrator_Module.Update_Init_Step!(integrator)
    ξ = Get_ξ(integrator)
    semi_implicit.wave_matrix .= Build_Wave_Matrices(semi_implicit.num_wavenumbers, integrator.laplacian_eigen, semi_implicit.div_mat, ξ)
end



"""
    Linear_Press_Gradient_δps(
        vert_coord, atmo_data,
        ps_ref, Δp_ref,
        lnp_half_ref, lnp_full_ref,
        t_ref
    )

Computes the reference state coefficient vector representing the sensitivity 
of the pressure gradient term to surface pressure variations (∂(RT ∇lnp)/∂pₛ).

### Parameters
    - vert_coord: Vertical grid configuration (providing bₖ coefficients).
    - atmo_data: Physical constants (gas constant R).

    - ps_ref: Reference surface pressure.
    - Δp_ref: Reference pressure thickness of layers.
    
    - lnp_half_ref: Reference log-pressure at interfaces.
    - lnp_full_ref: Reference log-pressure at layer centers.
    
    - t_ref: Reference temperature profile.

### Returns
    - R_T_δlnp: A component of the linearized gravity wave operator. 

"""
function Linear_Press_Gradient_δps(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    ps_ref::Float64, Δp_ref::Array{Float64,1},
    lnp_half_ref::Array{Float64,1}, lnp_full_ref::Array{Float64,1},
    t_ref::Array{Float64,1}
)
    
    nd, bk = vert_coord.nd, vert_coord.bk
    rdgas  = atmo_data.rdgas
    
    vert_difference_option = vert_coord.vert_difference_option
    
    R_T_δlnp = zeros(Float64, nd)

    if (vert_difference_option == "simmons_and_burridge")
        for k = 1:nd
            Δlnp_p = lnp_half_ref[k+1] - lnp_full_ref[k]
            Δlnp_m = lnp_full_ref[k]   - lnp_half_ref[k]
            R_T_δlnp[k] = rdgas * t_ref[k] * (bk[k+1] * Δlnp_p + bk[k] * Δlnp_m) / Δp_ref[k]
        end
    end

    return R_T_δlnp
end



"""
    Linear_Geopot_δt!(
        vert_coord, atmo_data,
        lnp_half_ref, lnp_full_ref,
        spe_δt, t_ref,
        spe_δgeopot_half, spe_δgeopot
    )

Integrates the hydrostatic equation to compute geopotential perturbations resulting 
from temperature perturbations (δΦ = γ δT).

### Parameters
    - vert_coord: Vertical grid configuration.
    - atmo_data: Physical constants (specifically gas constant R).
    
    - lnp_half_ref: Reference log-pressure at interfaces.
    - lnp_full_ref: Reference log-pressure at layer centers.
    
    - spe_δt: Spectral temperature perturbation input [nf, ns, nd].
    - t_ref: Reference temperature profile (passed for consistency, unused in linear term).
    
    - spe_δgeopot_half: Output geopotential perturbation at interfaces [nf, ns, nd+1].
    - spe_δgeopot: Output geopotential perturbation at layer centers [nf, ns, nd].

### Returns
    thing

### Modified
    - spe_δgeopot_half
    - spe_δgeopot: γ matrix in the semi-implicit system struct.

"""
function Linear_Geopot_δt!(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    lnp_half_ref::Array{Float64,1}, lnp_full_ref::Array{Float64,1},
    spe_δt::Array{ComplexF64,3}, t_ref::Array{Float64,1},
    spe_δgeopot_half::Array{ComplexF64,3}, spe_δgeopot::Array{ComplexF64,3}
)

    # compute δΦ = δΦ(t, lnp_half(ps), lnp_full(ps))
    #            = ∂Φ∂t(t_ref, lnp_half(ps)_ref, lnp_full(ps)_ref)δt 
    #            + [∂Φ∂lnp_half]_ref δlnp_half + [∂Φ∂lnp_full]_ref δlnp_full
    #            = ∂Φ∂t δt + ∂Φ∂ps ps_ref δlnps 
    # This function compute ∂Φ∂t δt

    ns, nf, nd = size(spe_δt)
    rdgas      = atmo_data.rdgas

    spe_δgeopot_half[:, :, nd+1] .= 0.0

    for k = nd:-1:2
        #todo optimize
        #Φ_{k-1/2} = Φ_{k+1/2} + RT_k(ln p_{k+1/2} - ln p_{k-1})
        spe_δgeopot_half[:, :, k] .= spe_δgeopot_half[:, :, k+1] + rdgas * (spe_δt[:, :, k] * (lnp_half_ref[k+1] - lnp_half_ref[k]))
    end

    for k = 1:nd
        #todo optimize
        #Φ_{k} = Φ_{k+1/2} + RT_k(ln p_{k+1/2} - ln p_{k})
        spe_δgeopot[:, :, k] .= spe_δgeopot_half[:, :, k+1] + rdgas * (spe_δt[:, :, k] * (lnp_half_ref[k+1] - lnp_full_ref[k]))
    end
end



"""
    Linear_Geopot_δps(
        vert_coord, atmo_data,
        δlnp_half_ref, δlnp_full_ref,
        lnp_half_ref, lnp_full_ref,
        t_ref
    )

Computes the reference state coefficient vector (H₂) representing the sensitivity 
of the geopotential height to surface pressure variations (∂Φ/∂pₛ).

### Parameters
    - vert_coord: Vertical grid configuration.
    - atmo_data: Physical constants (gas constant R).

    - δlnp_half_ref: Derivative of interface log-pressure w.r.t surface pressure (∂lnp_k+½/∂pₛ).
    - δlnp_full_ref: Derivative of layer log-pressure w.r.t surface pressure (∂lnp_k/∂pₛ).
    
    - lnp_half_ref: Reference log-pressure at interfaces
    - lnp_full_ref: Reference log-pressure at layer centers
    
    - t_ref: Reference temperature profile.

### Returns
    - δgeopot_ref: A vector of length [nd] (H₂) where the k-th element is the accumulated 
    geopotential response to a surface pressure perturbation.

"""
function Linear_Geopot_δps(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    δlnp_half_ref::Array{Float64,1}, δlnp_full_ref::Array{Float64,1},
    lnp_half_ref::Array{Float64,1}, lnp_full_ref::Array{Float64,1},
    t_ref::Array{Float64,1}
)

    # compute δΦ = δΦ(t, lnp_half(ps), lnp_full(ps))
    #            = ∂Φ∂t(t_ref, lnp_half(ps)_ref, lnp_full(ps)_ref)δt 
    #            + [∂Φ∂lnp_half]_ref δlnp_half + [∂Φ∂lnp_full]_ref δlnp_full
    #            = ∂Φ∂t δt + ∂Φ∂ps ps_ref δlnps 
    #
    # This function compute ∂Φ∂ps_ref

    nd = vert_coord.nd
    rdgas = atmo_data.rdgas
    δgeopot_ref = zeros(Float64, nd)
    δgeopot_half_ref = zeros(Float64, nd + 1)

    for k = nd:-1:2
        δgeopot_half_ref[k] = δgeopot_half_ref[k+1] + rdgas * (t_ref[k] * (δlnp_half_ref[k+1] - δlnp_half_ref[k]))
    end

    for k = 1:nd
        δgeopot_ref[k] = δgeopot_half_ref[k+1] + rdgas * (t_ref[k] * (δlnp_half_ref[k+1] - δlnp_full_ref[k]))
    end

    return δgeopot_ref
end


"""
    Linear_Ps_T_δdiv!(
        vert_coord, atmo_data, 
        spe_δdiv,
        ps_ref, Δp_ref,
        lnp_half_ref, lnp_full_ref,
        t_ref,
        spe_M_half, spe_Mdt_half,
        spe_δps, spe_δt
    )

Computes the linearized tendency of surface pressure (spe_δps) and temperature (spe_δt) 
resulting from divergence perturbations (spe_δdiv). This corresponds to calculating 
the terms -ν⋅δdiv and -τ⋅δdiv in the semi-implicit formulation.

### Parameters
    - vert_coord: Vertical grid configuration (Δa, Δb coefficients).
    - atmo_data: Physical constants (κ).

    - spe_δdiv: Input spectral divergence perturbation [nf, ns, nd].
    
    - ps_ref: Reference surface pressure.

    - Δp_ref: Reference pressure thickness of layers.
    
    - lnp_half_ref: Reference log-pressure at interfaces.
    - lnp_full_ref: Reference log-pressure at layer centers.
    
    - t_ref: Reference temperature profile.
    
    - spe_M_half: Workspace for vertical mass flux [nf, ns, nd+1].
    - spe_Mdt_half: Workspace for vertical advection term [nf, ns, nd+1].
    
    - spe_δps: Output surface pressure tendency [nf, ns, 1].
    - spe_δt: Output temperature tendency [nf, ns, nd].

### Returns
    - nothing

### Modified
    - spe_δps
    - spe_δt
    - spe_M_half
    - spe_Mdt_half

"""
function Linear_Ps_T_δdiv!(
    vert_coord::Vert_Coordinate, atmo_data::Atmo_Data, 
    spe_δdiv::Array{ComplexF64,3},
    ps_ref::Float64, Δp_ref::Array{Float64,1},
    lnp_half_ref::Array{Float64,1}, lnp_full_ref::Array{Float64,1},
    t_ref::Array{Float64,1},
    spe_M_half::Array{ComplexF64,3}, spe_Mdt_half::Array{ComplexF64,3},
    spe_δps::Array{ComplexF64,3}, spe_δt::Array{ComplexF64,3}
)
    # For temperature 
    # δ(-dσ ∂T∂σ + κTw/p) = -tau δdiv
    #                     = δ(-dσ∂T∂σ + κTw/p)
    #                     = δ(-dσ∂T∂σ) + [κT]_ref δ(w/p) (D_r depends on div)
    #                    
    # dmeam = D_k = div_k Δp_k, dmean_tot = ∑_{r=1}^{k-1} div_r Δp_r =  ∑_{r=1}^{k-1} Dr
    # w_k/p_k = dlnp/dt = ∂lnp/∂t + dσ ∂lnp/∂σ + v∇lnp
    #         = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k + v∇lnp
    #         = -[(∑_{r=1}^{k-1} Dr)(lnp_k+1/2 - lnp_k-1/2) + D_k(lnp_k+1/2 - lnp_k)]/Δp_k (∇lnp = 0)
    #
    # ∂ps/∂t = -∑ div_r Δp_r = -dmean_tot
    # M_{k+1/2} = -∑_{r=1}^k ∇(vrΔp_r) - B_{k+1/2}∂ps/∂t
    # The vertical discretization is 
    # [dσ∂T/∂σ]_k = 0.5(M_{k+1/2}(T_k+1 - T_k) + M_{k-1/2}(T_k - T_k-1))/Δp_k
    #
    # For surface pressure
    # δ(-∑∇(vk Δpk)) = δ(-∑div_k Δpk ) = -∑Δpk δdiv_k  = -nu δdiv
    #
    # div = [0,0, ..1,...0], spe_δps -> -nu, spe_δt ->  -tau[:,k]

    kappa      = atmo_data.kappa
    nf, ns, nd = size(spe_δdiv)
    vert_difference_option = vert_coord.vert_difference_option
    Δak, Δbk, bk = vert_coord.Δak, vert_coord.Δbk, vert_coord.bk

    dmean_tot = zeros(ComplexF64, nf, ns)
    dmean     = zeros(ComplexF64, nf, ns)

    if (vert_difference_option == "simmons_and_burridge")
        for k = 1:nd
            # Integrates the divergence column to find the total mass tendency
            @assert(Δak[k] + Δbk[k] * ps_ref ≈ Δp_ref[k])

            Δlnp_p = lnp_half_ref[k+1] - lnp_full_ref[k]
            Δlnp   = lnp_half_ref[k+1] - lnp_half_ref[k]
            
            # dmean = ∇ (vk Δp_k) = ∇vk Δp_k = Dk
            dmean .= spe_δdiv[:, :, k] * Δp_ref[k]

            spe_δt[:, :, k] .= -kappa * t_ref[k] * (dmean_tot * Δlnp + dmean * Δlnp_p) / Δp_ref[k]
            
            # dmean_tot = ∑_r=1^k Dr
            dmean_tot .+= dmean
            spe_M_half[:, :, k+1] .= -dmean_tot
        end
    end

    spe_δps[:, :, 1] .= -dmean_tot

    for k = 1:nd-1
        spe_M_half[:, :, k+1] .+= dmean_tot * bk[k+1]
    end

    spe_M_half[:, :, 1] .= 0.0
    spe_M_half[:, :, nd+1] .= 0.0

    # approximate the vertical advection term
    for k = 2:nd
        spe_Mdt_half[:, :, k] .= spe_M_half[:, :, k] * (t_ref[k] - t_ref[k-1])
    end

    spe_Mdt_half[:, :, 1] .= 0.0
    spe_Mdt_half[:, :, nd+1] .= 0.0
    for k = 1:nd
        spe_δt[:, :, k] .-= 0.5 * (spe_Mdt_half[:, :, k+1] + spe_Mdt_half[:, :, k]) / Δp_ref[k]
    end

end



"""
    Adjust_δlnps_δt_δdiv!(
        semi_implicit, vert_coord, atmo_data,
        spe_div_c, spe_div_p,
        spe_lnps_c, spe_lnps_p,
        spe_t_c, spe_t_p,
        spe_δdiv, spe_δlnps, spe_δt
    )

Updates the spectral tendencies for divergence, surface pressure, and temperature 
by adding the semi-implicit corrections. This step essentially prepares the 
RHS forcing for the implicit gravity wave solver.

### Parameters
    - semi_implicit: Solver structure containing reference arrays and temporary storage.
    - vert_coord: Vertical grid configuration.
    - atmo_data: Physical constants.

    - spe_div_c, spe_div_p: Spectral divergence at current (t) and previous (t-1) steps.

    - spe_lnps_c, spe_lnps_p: Spectral log-surface pressure at current and previous steps.

    - spe_t_c, spe_t_p: Spectral temperature at current and previous steps.

    - spe_δdiv: Input/Output spectral divergence tendency. Modified in-place.
    - spe_δlnps: Input/Output spectral log-pressure tendency. Modified in-place.
    - spe_δt: Input/Output spectral temperature tendency. Modified in-place.

### Returns
    - nothing

### Modified
    - spe_δdiv
    - spe_δlnps
    - spe_δt

See Implicit_Correction!
Adjust_δlnps_δt_δdiv!:  [f^p(i) , f^t(i),  f^d(i)] -> [Δlnps, Δt, Δdiv] := 
[f^p(i)+δlnps_temp , f^t(i)+δt_temp,  f^d(i) + I^d(T(i-1) - T(i) + ξ(f^t(i) + δt_temp), lnps(i-1) - lnps(i) + ξ(f^p(i) + δlnps_temp))] = 
[f^p(i)+δlnps_temp , f^t(i)+δt_temp,  f^d(i) + I^d(T(i-1) - T(i) + ξΔt, lnps(i-1) - lnps(i) + ξΔlnps)]

δlnps_temp = I^p(div(i-1) - div(i))
δt_temp = I^t(div(i-1) - div(i))

I^d(T(i-1) - T(i) + ξΔt, lnps(i-1) - lnps(i) + ξΔlnps)
= -∇^2(γ(T(i-1) - T(i) + ξΔt) + (H1+H2) ps_ref (lnps(i-1) - lnps(i) + ξΔlnps)))

"""
function Adjust_δlnps_δt_δdiv!(
    semi_implicit::Semi_Implicit_Solver, vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    spe_div_c::Array{ComplexF64,3}, spe_div_p::Array{ComplexF64,3},
    spe_lnps_c::Array{ComplexF64,3}, spe_lnps_p::Array{ComplexF64,3},
    spe_t_c::Array{ComplexF64,3}, spe_t_p::Array{ComplexF64,3},
    spe_δdiv::Array{ComplexF64,3}, spe_δlnps::Array{ComplexF64,3}, spe_δt::Array{ComplexF64,3}
)

    t_ref, ps_ref, Δp_ref, lnp_half_ref, lnp_full_ref = semi_implicit.t_ref, semi_implicit.ps_ref, semi_implicit.Δp_ref, semi_implicit.lnp_half_ref, semi_implicit.lnp_full_ref
    laplacian_eigen, h = semi_implicit.integrator.laplacian_eigen, semi_implicit.h

    spe_δdiv_temp, spe_δps_temp, spe_δt_temp, spe_δgeopot_half_temp, spe_δgeopot_temp = semi_implicit.spe_δdiv_temp,
    semi_implicit.spe_δps_temp, semi_implicit.spe_δt_temp, semi_implicit.spe_δgeopot_half_temp, semi_implicit.spe_δgeopot_temp

    # Computes the difference in state: Dᵖ - Dᶜ.
    spe_δdiv_temp .= spe_div_p - spe_div_c

    # spe_M_half::Array{ComplexF64,3}, spe_Mdt_half::Array{ComplexF64,3},
    spe_M_half_temp, spe_Mdt_half_temp = semi_implicit.spe_M_half_temp, semi_implicit.spe_Mdt_half_temp

    # Evaluates the linear response of temperature and pressure to this divergence difference using Linear_Ps_T_δdiv!.
    # Compute linearized δt, with constant surface pressure 
    Linear_Ps_T_δdiv!(vert_coord, atmo_data, spe_δdiv_temp, ps_ref, Δp_ref, lnp_half_ref, lnp_full_ref, t_ref,
        spe_M_half_temp, spe_Mdt_half_temp, spe_δps_temp, spe_δt_temp)

    spe_δt .+= spe_δt_temp
    spe_δlnps .+= spe_δps_temp / ps_ref

    # use as a memory container
    spe_δlnps_temp = semi_implicit.spe_δps_temp

    #todo
    ξ = Get_ξ(semi_implicit.integrator)

    spe_δt_temp .= spe_t_p - spe_t_c + ξ * spe_δt
    spe_δlnps_temp .= spe_lnps_p - spe_lnps_c + ξ * spe_δlnps

    # Updates the divergence tendency (spe_δdiv) by 
    # subtracting the Laplacian of the linear geopotential and pressure gradient terms derived from δT* and δlnpₛ*.
    # δD_final = δD_explicit - ∇²(Φ(δT*) + RT₀ ∇lnpₛ*)
    Linear_Geopot_δt!(vert_coord, atmo_data, lnp_half_ref, lnp_full_ref, spe_δt_temp, t_ref, spe_δgeopot_half_temp, spe_δgeopot_temp)

    nd = vert_coord.nd

    for k = 1:nd
        spe_δdiv[:, :, k] .-= laplacian_eigen .* (spe_δgeopot_temp[:, :, k] .+ h[k] * ps_ref * spe_δlnps_temp[:, :, 1])
    end

end



"""
    Implicit_Correction!(
            semi_implicit, vert_coord, atmo_data,
            spe_div_c, spe_div_p,
            spe_lnps_c, spe_lnps_p,
            spe_t_c, spe_t_p,
            spe_δdiv, spe_δlnps, spe_δt
        )

Solves the semi-implicit system in spectral space to filter fast-moving gravity waves. 
This function acts as the core solver, inverting the Helmholtz equation for divergence 
and back-substituting the results to update temperature and surface pressure tendencies.

### Parameters
    - semi_implicit: Solver structure containing the pre-computed inverse matrices.
    - vert_coord: Vertical grid configuration.
    - atmo_data: Physical constants.

    - spe_div_c, spe_div_p: Spectral divergence at current (t) and previous (t-1) steps.
    
    - spe_lnps_c, spe_lnps_p: Spectral log-surface pressure at current and previous steps.
    
    - spe_t_c, spe_t_p: Spectral temperature at current and previous steps.
    
    - spe_δdiv: Input/Output spectral divergence tendency.
    - spe_δlnps: Input/Output spectral log-pressure tendency.
    - spe_δt: Input/Output spectral temperature tendency.

### Returns
    - nothing

### Modified
    - spe_δdiv
    - spe_δlnps
    - spe_δt

The governing equations are
∂div/∂t  = ∇ × (A, B) - ∇^2E := f^d                    
∂lnps/∂t = (-∑_k div_k Δp_k + v_k ∇ Δp_k)/ps := f^p    
∂T/∂t    = -(u,v)∇T - dσ∂T∂σ + κTw/p + J := f^t           
Φ        = f^Φ                                               

implicit part: -∇^2Φ - ∇(RT∇lnp) ≈ I^d = -∇^2(γT + H2 ps_ref lnps) - ∇^2 H1 ps_ref lnps, here RT∇lnp ≈  H1 ps_ref ∇lnps
implicit part:  f^p              ≈ I^p = -ν div / ps_ref
implicit part:  - dσ∂T∂σ + κTw/p ≈ I^t = -τ div  
implicit part:  f^Φ              ≈ I^Φ = γT + H2 ps_ref lnps 
ν is a row vectors, H1, H2 are colume vectors, τ, γ and are matrices

We have 
δdiv  = f^d - I^d + I^d := E^d + I^d
δlnps = f^p - I^p + I^p := E^p + I^p
δT    = f^t - I^t + I^t := E^t + I^t

Leapfrog scheme  
δlnps = f^p(i) - I^p(i) + I^p(α div(i+1) + (1-α)div(i-1))
      = f^p(i) + I^p(div(i-1) - div(i)) + ξ I^p(δdiv) 
     := f^p(i) + δlnps_temp + ξ I^p(δdiv) 
      
δt    = f^t(i) - I^t(i) + I^t(α div(i+1) + (1-α)div(i-1))
      = f^t(i) + I^t(div(i-1) - div(i)) + ξ I^t(δdiv)
     := f^t(i) + δt_temp + ξ I^t(δdiv)

δdiv  = f^d(i) - I^d(i) + I^d(α T(i+1)+(1-α)T(i-1), α lnps(i+1) + (1-α)lnps(i-1))
      = f^d(i) + I^d(T(i-1) - T(i), lnps(i-1) - lnps(i)) + ξ I^d(δT, δlnps)
      = f^d(i) + I^d(T(i-1) - T(i), lnps(i-1) - lnps(i)) + ξ I^d(f^t(i) + δt_temp + I^t(ξ δdiv), f^p(i) + δlnps_temp + ξ I^p(δdiv))
      = f^d(i) + I^d(T(i-1) - T(i) + ξ(f^t(i) + δt_temp), lnps(i-1) - lnps(i) + ξ(f^p(i) + δlnps_temp)) + ξ I^d( ξ I^t( δdiv),  ξ I^p(δdiv))

For the first step, Backward-Euler scheme
set Var(i-1) = Var(i), and ξ = αΔt


Adjust_δlnps_δt_δdiv!:  [f^p(i) , f^t(i),  f^d(i)] -> [Δlnps, Δt, Δdiv] := [f^p(i)+δlnps_temp , f^t(i)+δt_temp,  f^d(i) + I^d(T(i-1) - T(i) + ξ(f^t(i) + δt_temp), lnps(i-1) - lnps(i) + ξ(f^p(i) + δlnps_temp))]
Then we have 

δlnps = Δlnps + ξ I^p(δdiv) 
δt    = Δt + ξ I^t(δdiv)
δdiv  = Δdiv + ξ^2 I^d( I^t( δdiv),  I^p(δdiv))

Solve for δdiv, δt, δlnps sequentially
    δdiv = Δdiv - ξ^2 ∇^2(γ I^t( δdiv) + (H2+H1) ps_ref I^p(δdiv))
         = Δdiv - ξ^2 ∇^2(γ (-τ) δdiv + (H2+H1) ps_ref (-ν) δdiv / ps_ref)
         = Δdiv - ξ^2 ∇^2(γ (-τ) δdiv + (H2+H1)  (-ν) δdiv)
    (1 - ξ^2 ∇^2(γ τ + (H2+H1)ν)δdiv = Δdiv
  
  
    wave_matrix[:,:,s] = (I + ξ^2 (γ τ + (H2+H1)ν) s(s+1)/r^2)^{-1}

"""
function Implicit_Correction!(
    semi_implicit::Semi_Implicit_Solver, vert_coord::Vert_Coordinate, atmo_data::Atmo_Data,
    spe_div_c::Array{ComplexF64,3}, spe_div_p::Array{ComplexF64,3},
    spe_lnps_c::Array{ComplexF64,3}, spe_lnps_p::Array{ComplexF64,3},
    spe_t_c::Array{ComplexF64,3}, spe_t_p::Array{ComplexF64,3},
    spe_δdiv::Array{ComplexF64,3}, spe_δlnps::Array{ComplexF64,3}, spe_δt::Array{ComplexF64,3}
)

    ξ = Get_ξ(semi_implicit.integrator)

    # Calls Adjust_δlnps_δt_δdiv! to prepare the forcing terms for the Helmholtz equation.
    # This adds the gravity wave contributions from the previous time steps to the explicit tendencies.
    Adjust_δlnps_δt_δdiv!(semi_implicit, vert_coord, atmo_data,
        spe_div_c, spe_div_p,
        spe_lnps_c, spe_lnps_p,
        spe_t_c, spe_t_p,
        spe_δdiv, spe_δlnps, spe_δt)


    num_wavenumbers = semi_implicit.num_wavenumbers
    wave_numbers = semi_implicit.wave_numbers
    wave_matrix = semi_implicit.wave_matrix

    # Helmholtz Inversion (Divergence Solve)
    # Solves for the implicit divergence tendency (δD) for each total wavenumber L:
    # δDₗ = (I - ξ² B ∇²)⁻¹ ⋅ RHSₗ
    # where B is the vertical structure matrix (div_mat) and the inverse operator is stored in wave_matrix.
    nf, ns, nd = size(spe_δdiv)
    for m = 0:nf-1
        for n = m:ns-1
            L = wave_numbers[m+1, n+1]
            @assert(L == n)
            # does not need the last spherical mode
            if (L <= num_wavenumbers)
                spe_δdiv[m+1, n+1, :] .= wave_matrix[:, :, L+1] * spe_δdiv[m+1, n+1, :]
            end
        end
    end

    t_ref, ps_ref, Δp_ref, lnp_half_ref, lnp_full_ref = semi_implicit.t_ref, semi_implicit.ps_ref, semi_implicit.Δp_ref, semi_implicit.lnp_half_ref, semi_implicit.lnp_full_ref

    spe_δps_temp, spe_δt_temp = semi_implicit.spe_δps_temp, semi_implicit.spe_δt_temp
    spe_M_half_temp, spe_Mdt_half_temp = semi_implicit.spe_M_half_temp, semi_implicit.spe_Mdt_half_temp

    # Back-Substitution
    # Updates the temperature and surface pressure tendencies using the newly solved 
    # divergence tendency:
    # δT = δT* + ξ ⋅ (-τ ⋅ δD)
    # δlnpₛ = δlnpₛ* + ξ ⋅ (-ν ⋅ δD)
    # This ensures the final tendencies are consistent with the semi-implicit formulation.
    Linear_Ps_T_δdiv!(vert_coord, atmo_data, spe_δdiv,
        ps_ref, Δp_ref, lnp_half_ref, lnp_full_ref, t_ref,
        spe_M_half_temp, spe_Mdt_half_temp,
        spe_δps_temp, spe_δt_temp)

    spe_δt .+= ξ * spe_δt_temp
    spe_δlnps .+= ξ / ps_ref * spe_δps_temp

end

end