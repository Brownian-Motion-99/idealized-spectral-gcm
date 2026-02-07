module Time_Integrator_Module

export Filtered_Leapfrog, Update_Init_Step!, Get_Δt, Get_ξ, Compute_Spectral_Damping!, Filtered_Leapfrog!


mutable struct Filtered_Leapfrog
    robert_coef::Float64

    damping_order::Int64
    damping_coef::Float64
    damping::Array{Float64,2}
    laplacian_eigen::Array{Float64,2}

    implicit_coef::Float64

    Δt::Int64
    init_step::Bool
    time::Int64
    start_time::Int64
    end_time::Int64

end



"""
    Filtered_Leapfrog(
        robert_coef, damping_order, damping_coef, 
        eigen,
        implicit_coef,
        Δt, init_step, start_time, end_time
    )

Initializes the time-integration state and pre-calculates spectral damping coefficients 
based on the Laplacian eigenvalues. This constructor prepares the operator for a 
Robert-Asselin-Williams filtered leapfrog scheme combined with scale-selective 
hyper-diffusion.

### Parameters
    - robert_coef: The Robert-Asselin filter coefficient (ν_RA) used to damp the computational mode.
    - damping_order: The order of the diffusion operator (2n). Must be an even integer.
    - damping_coef: The hyper-diffusion strength coefficient (K).

    - eigen: Eigenvalues (λ) of the Laplacian operator corresponding to the spectral/grid resolution.
    
    - implicit_coef: Off-centering parameter for semi-implicit discretization (α).
    
    - Δt: The fundamental time step increment in seconds.
    - init_step: Boolean flag indicating if the integration is at the first step.
    - start_time: Initial simulation time in seconds.
    - end_time: Final simulation time in seconds.

### Returns
    - Filtered_Leapfrog: An initialized mutable structure containing the integration state and damping matrix.

"""
function Filtered_Leapfrog(
    robert_coef::Float64, damping_order::Int64, damping_coef::Float64, 
    eigen::Array{Float64, 2},
    implicit_coef::Float64,
    Δt::Int64, init_step::Bool, start_time::Int64, end_time::Int64
)

    @assert(damping_order%2 == 0)
    
    num_fourier, num_spherical = size(eigen) .- 1 
    
    # resolution_independent damping
    # damping = damping_coef*(eigen).^damping_order
    # resolution_dependent damping

    # The spectral damping is normalized relative to the maximum eigenvalue: D = K ⋅ (λ / λ_max)^n.
    damping = damping_coef * (eigen / eigen[1,num_spherical]).^damping_order

    laplacian_eigen = eigen
    time = start_time
    Filtered_Leapfrog(
        robert_coef, damping_order, damping_coef, damping, 
        laplacian_eigen, 
        implicit_coef, 
        Δt, init_step, time, start_time, end_time
    )
end



function Update_Init_Step!(integrator::Filtered_Leapfrog)
    # @assert(integrator.init_step)
    integrator.init_step = false
end



function Get_Δt(integrator::Filtered_Leapfrog)
    init_step = integrator.init_step
    Δt = integrator.Δt
    return (init_step ? Δt : 2*Δt) 
end



function Get_ξ(integrator::Filtered_Leapfrog)
    init_step = integrator.init_step
    Δt, implicit_coef = integrator.Δt, integrator.implicit_coef
    return (init_step ? Δt*implicit_coef : 2*Δt*implicit_coef) 
end



"""
    Compute_Spectral_Damping!(integrator, Qc, Qp, δQ)

Updates the spectral tendency (δQ) by applying a scale-selective hyper-diffusion operator. 
The implementation uses an implicit formulation for numerical stability, particularly 
to handle the stiffness associated with high-order Laplacian operators.

### Parameters
    - integrator: The Filtered_Leapfrog structure containing damping coefficients and time-step state.
    - Qc: Current time-level spectral coefficients (Qⁱ).
    - Qp: Previous time-level spectral coefficients (Qⁱ⁻¹).
    - δQ: The spectral tendency (dQ/dt) to be modified by the damping term.

### Returns
    - nothing

### Modified
    - δQ

"""
function Compute_Spectral_Damping!(
    integrator::Filtered_Leapfrog,
    Qc::AbstractArray{ComplexF64,3}, Qp::AbstractArray{ComplexF64,3}, 
    δQ::AbstractArray{ComplexF64,3}
)

    # (Q^{i+1} - Q^{i-1}) / 2Δt = dQ - ν(-1)^n ∇^2n Q^{i+1}
    #                           = dQ - ν(-1)^n ∇^2n(Q^{i+1} -  Q^{i-1}) - ν(-1)^n ∇^2n Q^{i-1}

    # (Q^{i+1} - Q^{i-1}) (1/2Δt + ν |σ|^n) = dQ - ν|σ|^n  Q^{i-1}

    # (Q^{i+1} - Q^{i-1})/2Δt = (dQ - ν|σ|^n Q^{i-1}) / (1 +  ν 2Δt |σ|^n)

    init_step = integrator.init_step
    damping = integrator.damping
    Δt = integrator.Δt
      
    if (init_step) 
        δQ .= (δQ - Qc .* damping) ./ (1.0 .+ Δt*damping)
    else
        δQ .= (δQ - Qp .* damping) ./ (1.0 .+ 2Δt*damping)
    end
end



"""
    Filtered_Leapfrog!(integrator, δQ, Qp, Qc, Qn)

Performs the temporal update of spectral coefficients using a three-level leapfrog 
scheme integrated with a Robert-Asselin filter (RAW filter). The function handles 
both the initial forward-Euler step and the subsequent filtered leapfrog steps.

### Parameters
    - integrator: The Filtered_Leapfrog structure containing coefficients and time state.
    - δQ: The spectral tendency (dQ/dt), typically pre-processed by damping.
    - Qp: Spectral coefficients at the previous time level (Qⁱ⁻¹).
    - Qc: Spectral coefficients at the current time level (Qⁱ).
    - Qn: Spectral coefficients at the next time level (Qⁱ⁺¹).

### Returns
    - nothing

### Modified
    - Qc
    - Qn

"""
function Filtered_Leapfrog!(
    integrator::Filtered_Leapfrog,
    δQ::AbstractArray{ComplexF64,3},
    Qp::AbstractArray{ComplexF64,3}, Qc::AbstractArray{ComplexF64,3}, Qn::AbstractArray{ComplexF64,3}
)

    init_step = integrator.init_step
    robert_coef = integrator.robert_coef
    Δt = integrator.Δt
    
    if (init_step) 
        Qn  .= Qc + Δt * δQ
        Qc .+= robert_coef * (-1.0 * Qc + Qn)    
    else
        Qc .+= robert_coef * (Qp - 2 * Qc)
        Qn  .= Qp + 2 * Δt * δQ
        Qc .+= robert_coef * Qn
    end

end

end