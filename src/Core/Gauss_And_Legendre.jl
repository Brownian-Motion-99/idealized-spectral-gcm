module Gauss_And_Legendre_Module

using LinearAlgebra

export Compute_Legendre, Compute_Gaussian



"""
    Compute_Legendre(num_fourier, num_spherical, sinθ, nθ)

Computes the Associated Legendre Polynomials and their first derivatives with respect to
sine of latitude (μ) using stable three-term recurrence relations.

Mathematical Formulation:
1. Normalization: P_{0,0} = 1
2. Diagonal Recurrence (Sectorial): P_{m,m} = √((2m+1)/2m) cosθ P_{m-1,m-1}
3. Sub-diagonal Recurrence: P_{m+1,m} = √(2m+3) sinθ P_{m,m}
4. Vertical Recurrence (General): ε_{l,m} P_{l,m} = sinθ P_{l-1,m} - ε_{l-1,m} P_{l-2,m}, where ε_{l,m} = √((l² - m²) / (4l² - 1))
5. Derivative (Meridional): (1-μ²) dP_{l,m}/dμ = -l ε_{l+1,m} P_{l+1,m} + (l+1) ε_{l,m} P_{l-1,m}

### Parameters
    - num_fourier: Maximum Fourier wavenumber
    - num_spherical: Maximum spherical wavenumber
    - sinθ: Sine values of Gaussian latitudes
    - nθ: Latitudinal grid size

### Returns
    - qnm (num_fourier+1, num_spherical+1, nθ): Associated Legendre Polynomials P_{l,m}(sinθ).
        Indexing: qnm[m+1, l+1, lat] maps to physical mode P_{l,m}.

    - dqnm (num_fourier+1, num_spherical+1, nθ): Meridional derivatives dP_{l,m}/dμ.
        Indexing: dqnm[m+1, l+1, lat] maps to physical derivative dP_{l,m}/dμ.

### Modified
    - nothing

"""
function Compute_Legendre(num_fourier, num_spherical, sinθ, nθ)

    # Associated Legendre polynomials and their derivatives
    qnm  = zeros(Float64, num_fourier+1, num_spherical+2, nθ)
    dqnm = zeros(Float64, num_fourier+1, num_spherical+1, nθ)
    
    cosθ = sqrt.(1 .- sinθ.^2)
    ε    = zeros(Float64, num_fourier+1, num_spherical+2)

    # The diagonal recurrence (l == m)
    qnm[1, 1, :] .= 1.0
    for m = 1:num_fourier
        qnm[m+1, m+1, :] = sqrt((2m+1)/(2m)) .* cosθ .* qnm[m, m, :]
    end
    
    # The semi-diagonal recurrence (l == m+1)
    for m = 1:num_fourier+1
        qnm[m, m+1, :] = sqrt(2*m+1) * sinθ .* qnm[m, m, :] 
    end
    
    # Normalization factors
    for m = 0:num_fourier
        for l = m:num_spherical+1
            ε[m+1, l+1] = sqrt((l^2 - m^2) ./ (4*l^2 - 1))
        end
    end

    # The main loop (l > m+1)
    for m = 0:num_fourier
        for l = m+2:num_spherical+1
            qnm[m+1, l+1, :] = (sinθ .* qnm[m+1, l, :] -  ε[m+1, l] * qnm[m+1, l-1, :]) / ε[m+1, l+1]
        end
    end

    # Derivatives
    for m = 0:num_fourier
        for l = m:num_spherical
            if l == m
                dqnm[m+1, l+1, :] = (-l * ε[m+1, l+2] * qnm[m+1, l+2, :]) ./ (cosθ.^2)
            else
                dqnm[m+1, l+1, :] = (-l * ε[m+1, l+2] * qnm[m+1, l+2, :] + (l+1) * ε[m+1, l+1] * qnm[m+1, l, :]) ./ (cosθ.^2)
            end
        end
    end

    return qnm[:, 1:num_spherical+1, :], dqnm
end
    


"""
    Compute_Gaussian(n)

Computes the Gaussian latitudes and quadrature weights for a spectral grid.

Mathematical Formulation:
1. Symmetry: The roots are symmetric about the equator (x=0). The function solves for the northern hemisphere roots and mirrors them.
2. Recurrence Relation: n Pn(x) = (2n-1) x P_{n-1}(x) - (n-1) P_{n-2}(x)
3. Derivative: (x² - 1) P'n(x) = n (x Pn(x) - P_{n-1}(x))
4. Newton-Raphson Update: x_{new} = x_{old} - Pn(x) / P'n(x)
5. Initial Guess (Asymptotic): x₀ = cos(π(i - 0.25) / (n + 0.5))
6. Quadrature Weights: wᵢ = 2 / ((1 - xᵢ²) [P'n(xᵢ)]²)

### Parameters
    - n: Latitudinal grid size (must be even)

### Returns
    - sinθ (n): Gaussian nodes (roots of Legendre polynomial), representing sine of latitude.
        Ordered from south (-z) to north (+z) if using the provided symmetry logic.

    - wts (n): Gaussian quadrature weights corresponding to the nodes in sinθ.

### Modified
    - nothing

"""
function Compute_Gaussian(n)
    
    itermax = 10000
    tol     = 1.0e-15

    sinθ = zeros(Float64, n)
    wts  = zeros(Float64, n)

    n_half = Int64(n/2)
    for i = 1:n_half
        dp = 0.0
        z = cos(pi*(i - 0.25)/(n + 0.5))
        
        for iter = 1:itermax
            p2 = 0.0
            p1 = 1.0
            
            for j = 1:n
                p3 = p2 # Pj-2
                p2 = p1 # Pj-1
                p1 = ((2.0*j - 1.0)*z*p2 - (j - 1.0)*p3)/j  #Pj
            end

            # P'_n
            dp = n * (z * p1 - p2) / (z * z - 1.0)
            z1 = z
            z  = z1 - p1 / dp
            if(abs(z - z1) <= tol)
                break;
            end
            if iter == itermax
                @error("Compute_Gaussian does not converge!")
            end
        end
        
        sinθ[i], sinθ[n-i+1] = -z, z
        wts[i] = wts[n-i+1]  = 2.0 / ((1.0 - z * z) * dp * dp)
    end

    return sinθ, wts
end

end