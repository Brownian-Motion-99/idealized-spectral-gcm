# Model Physics & Formulation

idealized-spectral-gcm solves the hydrostatic primitive equations on a sphere using the spectral transform method. The model separates the **Dynamical Core** (resolved flow) from **Physical Parameterizations** (sub-grid processes like large scale condensation and boundary layer friction).

## Governing Equations (Dynamical Core)

The model integrates the dry or moist primitive equations using the **Vorticity-Divergence** formulation.

The prognostic variables are:
* Absolute Vorticity: $\eta = \zeta + f$
* Divergence: $D = \nabla \cdot \mathbf{v}$
* Temperature: $T$
* Log-Surface Pressure: $\Pi = \ln(p_s)$
* Specific Humidity: $q$ (optional)

### Momentum Equations
The model evolves the vertical component of absolute vorticity and horizontal divergence:

```math
\frac{\partial \eta}{\partial t} = -\frac{1}{a \cos^2\theta} \left[ \frac{\partial (U\eta)}{\partial \lambda} + \cos\theta \frac{\partial (V\eta)}{\partial \theta} \right] - \nabla \cdot (\mathbf{v} \dot{\sigma} \frac{\partial \mathbf{v}}{\partial \sigma})
```

```math
\frac{\partial D}{\partial t} = \nabla^2 (\Phi + E) - \nabla \cdot (\eta \mathbf{k} \times \mathbf{v}) - \nabla \cdot (\mathbf{v} \dot{\sigma} \frac{\partial \mathbf{v}}{\partial \sigma})
```

Where $E = \frac{u^2 + v^2}{2}$ is the kinetic energy per unit mass. $\Phi$ is the geopotential. $U = u \cos \theta, V = v \cos \theta$.

### Thermodynamic & Continuity Equations
```math
\frac{\partial T}{\partial t} = - \nabla \cdot (\mathbf{v} T) + T D + \dot{\sigma} \frac{\partial T}{\partial \sigma} + \frac{\kappa T \omega}{p} + \frac{Q}{c_p}
```

```math
\frac{\partial \Pi}{\partial t} = - \mathbf{v} \cdot \nabla \Pi - D - \frac{\partial \dot{\sigma}}{\partial \sigma}
```

---

## Semi-Implicit Time Integration

The model uses a **Semi-Implicit Leapfrog** scheme to maintain numerical stability with large time steps ($\Delta t$). Fast-moving gravity waves are treated implicitly, while advection and Coriolis terms are treated explicitly.

### Linearization
The state variables are linearized around a reference state of rest ($\bar{T}(\sigma)$, $\bar{p}_s$). The equations are split into linear terms (treated implicitly) and non-linear residuals ($N$):

```math
\frac{\partial D}{\partial t} + \nabla^2 (\Phi' + R \bar{T} \Pi') = N_D
```

```math
\frac{\partial T}{\partial t} + \tau D = N_T
```

```math
\frac{\partial \Pi}{\partial t} + \nu D = N_\Pi
```

Where $\tau$ and $\nu$ are vertical coupling vectors derived from the reference profile.

### The Helmholtz Equation
By eliminating $T$ and $\Pi$ in favor of $D$, the system reduces to a **Helmholtz equation** for the divergence tendency $\delta D$ in spectral space:

```math
(I - \xi^2 \mathbf{B} \nabla^2) \delta D = \text{RHS}
```

* **$\xi$:** Time step scaling factor ($\alpha \Delta t$).
* **$\mathbf{B}$:** A vertical structure matrix combining the hydrostatic and continuity effects (`div_mat` in the code).
* **$\nabla^2$:** Laplacian operator (diagonal in spectral space: $-n(n+1)/a^2$).

This system is solved using the pre-computed inverse matrices in `Semi_Implicit.jl`.

---

## Physical Parameterizations

The `physics_params` dictionary controls which sub-grid processes are active.

### Held-Suarez Forcing (`HS_Forcing.jl`)
Used for idealized dry dynamics benchmarks. It approximates radiative cooling and boundary layer friction.

* **Newtonian Relaxation (Temperature):**
    Relaxation towards a zonally symmetric equilibrium temperature $T_{eq}$:

$$\frac{\partial T}{\partial t} = -k_T (\phi, \sigma) (T - T_{eq})$$
    
* **Rayleigh Damping (Wind):**
    Linear damping of winds near the surface ($\sigma > \sigma_b$) to represent surface friction:
    
$$\frac{\partial \mathbf{v}}{\partial t} = -k_v (\sigma) \mathbf{v}$$

* **Frictional Heating**
    The kinetic energy loss due to Rayleigh Damping is returned to the system as sensible heat.

$$\frac{\partial T}{\partial t}_{fric} = -\frac{1}{c_p} (u \frac{\partial u}{\partial t}_{fric} + v \frac{\partial v}{\partial t}_{fric})$$

### Large-Scale Condensation (`Lscale_Cond.jl`)
A saturation adjustment diagnosed after convection. It uses the same mixed-phase
Smithsonian saturation-vapor-pressure equation and exact saturation-specific-humidity
conversion as Betts-Miller and the surface evaporation scheme. If specific
humidity $q$ exceeds $q_{sat}(T,p)$, the signed humidity increment is

```math
\Delta q_{LS} =
\frac{q_{sat}-q}
{1 + \frac{L_v}{c_p} \frac{\partial q_{sat}}{\partial T}} \leq 0.
```

The temperature increment is

```math
\Delta T_{LS} = -L\frac{L_v}{c_p}\Delta q_{LS} \geq 0,
```

where $L$ scales latent heating only. Humidity removal and precipitation do not
depend on $L$. For $L<1$, reduced heating is intentional, so exact moist-energy
closure and exact final saturation are not expected.

### Betts-Miller Convection (`Betts_Miller.jl`)

For a column with positive contiguous CAPE, the scheme constructs a moist
adiabatic parcel temperature profile and a reference humidity profile equal to
the configured relative humidity times the parcel saturation mixing ratio. It
returns the relaxation rates

```math
\dot T = \frac{T_{ref}-T}{\tau_{BM}}, \qquad
\dot q = \frac{q_{ref}-q}{\tau_{BM}}.
```

The column energy correction lengthens whichever of the temperature or
humidity adjustments supplies the larger precipitation-equivalent energy; it
does not cap either tendency using the model time step. Physics is subcycled at
the base model timestep, so the driver requires `Δt <= bm_tau`; each explicit
BM adjustment therefore cannot overshoot its reference profile.

Saturation vapor pressure follows the reference mixed-phase Smithsonian
thermodynamics: ice below 253.16 K, liquid above 273.16 K, and a linear blend
between them. The parcel calculation includes an initially saturated
surface-parcel adjustment, an RK2 saturated ascent using the arithmetic
midpoint pressure, and discrete contiguous LFC/LZB detection. Convection that
remains buoyant through the model top uses the top full level as its LZB.

When both moist schemes are enabled, Betts-Miller updates the private physics
state first and large-scale condensation diagnoses that updated state. The
complete physics chain is applied directly to the provisional post-dynamics
state in 600 s substeps: one substep during startup and two across a mature
leapfrog interval. Physics increments are not inserted into the leapfrog RHS.
The process order is BM, large-scale condensation, surface exchange, PBL
mixing, Rayleigh friction, Newtonian relaxation, and LRF. Diagnostic tendency
and flux fields are time averages across the substeps.

### Boundary Layer Mixing (`PBL.jl`)
Vertical turbulent diffusion of heat, moisture, and momentum.

* **Surface Fluxes:** Implicitly computed using bulk aerodynamic formulas:

```math
SH = \rho C_H |V_s| (T_{surf} - T_{air})
```
    
```math
LH = \rho C_E |V_s| (q_{surf} - q_{air})
```
    
* **Vertical Diffusion:** Solved implicitly using a tridiagonal solver (Thomas Algorithm):

```math
\frac{\partial \psi}{\partial t} = \frac{\partial}{\partial p} \left( g^2 \rho^2 K_E \frac{\partial \psi}{\partial p} \right)
```
