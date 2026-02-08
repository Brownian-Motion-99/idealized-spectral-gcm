# Numerics & Discretization

This section details the numerical methods used to discretize the primitive equations in space. The model uses a hybrid approach: a **Spectral Transform Method** for the horizontal directions and **Finite Differences** for the vertical direction.

## Horizontal Discretization: The Spectral Transform

The model represents the horizontal state of the atmosphere using truncated series of **Spherical Harmonics** $Y_n^m(\lambda, \mu)$. This allows for the exact calculation of horizontal derivatives on a sphere without pole problems.

### Spherical Harmonic Expansion
Any scalar field $\psi(\lambda, \mu)$ (where $\mu = \sin\theta$) is expanded as 

```math
\psi(\lambda, \mu) = \sum_{m=-M}^{M} \sum_{n=|m|}^{N(m)} \psi_n^m Y_n^m(\lambda, \mu)
```

The basis functions are defined as $Y_n^m(\lambda, \mu) = P_n^m(\mu) e^{im\lambda}$, where $P_n^m(\mu)$ are the **Associated Legendre Polynomials** of the first kind.

### The Gaussian Grid
To compute non-linear terms (like advection) efficiently, the model transforms variables from spectral space to a physical **Gaussian Grid**.

* **Longitude ($\lambda$):** Uniformly spaced points. Transforms use the Fast Fourier Transform (FFT).
* **Latitude ($\theta$):** Spaced according to the roots of the Legendre polynomial $P_J(\mu)$.

The quadrature weights $w_j$ satisfy the Gaussian quadrature rule, allowing exact integration of polynomials up to degree $2J - 1$ via the summation $\int_{-1}^{1} P(\mu) d\mu = \sum_{j=1}^{J} w_j P(\mu_j)$.

### Transform Cycle
The `Spectral_Spherical_Mesh` module manages the cycle between spaces:

* **Synthesis (Spectral to Grid):** Inverse Legendre Transform followed by Inverse FFT.
* **Non-linear Operations:** Computed locally on the grid (e.g., $u \cdot T$).
* **Analysis (Grid to Spectral):** Forward FFT followed by Forward Legendre Transform.

## Vertical Discretization: Finite Differences

The vertical structure uses a **Hybrid $\sigma-p$ Coordinate** system, discretized using the energy-conserving finite difference scheme of **Simmons and Burridge (1981)**.

### Vertical Coordinate Definition
The pressure at a model interface level $k+1/2$ is defined by two coefficients, $A_k$ and $B_k$, such that 

```math
p_{k+1/2} = A_{k+1/2} + B_{k+1/2} p_s
```

* **Top of Atmosphere:** $p_{1/2} = 0$ (Pure pressure).
* **Surface:** $p_{N+1/2} = p_s$ (Pure sigma).

### Grid Staggering (Lorenz Grid)
The model uses a standard Lorenz staggering to prevent decoupling of the pressure and velocity fields:

* **Full Levels ($k$):** $u, v, T, q, \omega$
* **Half Levels ($k+1/2$):** $\dot{\sigma}$ (vertical velocity), Fluxes.

### Hydrostatic Integration
Geopotential $\Phi$ is integrated upward from the surface geopotential $\Phi_s$:

```math
\Phi_k = \Phi_s + \sum_{j=k}^{N} R T_j \Delta \ln p_j
```

The term $\Delta \ln p$ is handled specifically by the Simmons-Burridge scheme to ensure that the discrete analytic relation $\nabla \Phi + RT \nabla \ln p$ balances exactly for the reference atmosphere, preventing spurious generation of angular momentum.

### Vertical Advection
Vertical advection of a scalar $\chi$ is computed in flux form to conserve mass: 

```math
(\dot{\sigma} \frac{\partial \chi}{\partial \sigma})_k \approx \frac{1}{\Delta p_k} [ (\dot{\sigma} \chi)_{k+1/2} - (\dot{\sigma} \chi)_{k-1/2} ]
```

The value of $\chi$ at the interface $k+1/2$ is interpolated from the full levels using a weighting scheme (linear or logarithmic) defined in `Vert_Coordinate.jl`.