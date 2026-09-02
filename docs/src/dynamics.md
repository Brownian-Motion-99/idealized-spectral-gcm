# Dynamical core

This page describes the resolved equations and numerical methods implemented in
`src/Dynamics` and `src/Core`. Parameterized tendencies are documented
separately in [Physical parameterizations](@ref).

The driver selects one of three equation sets with `model_type`:

| `model_type` | Prognostic system | Vertical levels |
|:---|:---|:---|
| `:PrimitiveEquation` | hydrostatic primitive equations | `nd` layers |
| `:Shallow_Water` | rotating shallow-water equations | one layer |
| `:Barotropic` | nondivergent barotropic vorticity equation | one layer |

All three systems use the same spherical-harmonic transforms, Gaussian grid,
leapfrog integrator, Robert--Asselin filter, and scale-selective spectral
damping. See [Notation](@ref) for the symbols and their code names.

## Hydrostatic primitive equations

### Prognostic and diagnostic variables

At each full model level the dry part of the core advances spectral
coefficients of relative vorticity $\zeta$, divergence $\delta$, and
temperature $T$. Log surface pressure $\Pi=\ln p_s$ is a two-dimensional
spectral prognostic. Specific humidity $q$ is optional and is advanced on the
Gaussian grid by a separate conservative finite-volume transport scheme.

The wind components $u$ and $v$, pressure $p$, geopotential $\Phi$, pressure
vertical velocity $\omega$, absolute vorticity $\eta=\zeta+f$, and interface
mass flux $M$ are diagnosed. This distinction matters when reading the source:
`spe_vor_*` stores **relative**, not absolute, vorticity, and `grid_w_full`
stores $\omega=Dp/Dt$ in Pa s$^{-1}$ rather than geometric vertical velocity.

For moist simulations the equation-of-state temperature is

```math
T_v = T\left[1 + \left(\frac{R_v}{R_d}-1\right)q\right].
```

$T_v$ is used consistently in hydrostatic integration, the pressure-gradient
term, pressure work, density estimates, and surface-layer height. Setting
`moisture_processes = true` therefore requires
`physics_params["use_virtual_temperature"] = true`.

### Vector-invariant momentum equations

The horizontal momentum equation is evaluated in vector-invariant form. Its
resolved, adiabatic part can be written schematically as

```math
\frac{\partial\mathbf{v}}{\partial t}
= -\eta\,\mathbf{k}\times\mathbf{v}
  -\nabla_h(\Phi+K)
  -R_dT_v\nabla_h\ln p
  -\mathcal{A}_v(\mathbf{v}),
```

where $K=(u^2+v^2)/2$ and $\mathcal{A}_v$ is the discrete vertical-advection
operator. The code first forms the gridpoint vector

```math
(F_u,F_v)=(\eta v,-\eta u)
```

together with the hybrid-coordinate pressure-gradient and vertical-advection
terms. `Vor_Div_From_Grid_UV!` then takes the rotational and divergent parts of
that vector to form the tendencies of $\zeta$ and $\delta$. The Laplacian of
$\Phi+K$ is added to the divergence equation spectrally, where the Laplacian is
diagonal.

This formulation avoids prognosing winds directly in spectral space and
retains a regular vorticity--divergence representation at the poles.

### Mass continuity and hybrid-coordinate mass flux

The interface pressures are

```math
p_{k+1/2}=A_{k+1/2}+B_{k+1/2}p_s,
\qquad
\Delta p_k=p_{k+1/2}-p_{k-1/2}.
```

Vertically integrating mass continuity gives the surface-pressure tendency

```math
\frac{\partial p_s}{\partial t}
=-\sum_{k=1}^{N}\nabla_h\!\cdot
  \left(\mathbf{v}_k\Delta p_k\right).
```

For each interface, the downward mass flux relative to the moving hybrid
surface is diagnosed as

```math
M_{k+1/2}
=-\sum_{r=1}^{k}\nabla_h\!\cdot
  \left(\mathbf{v}_r\Delta p_r\right)
-B_{k+1/2}\frac{\partial p_s}{\partial t},
```

with $M_{1/2}=M_{N+1/2}=0$. `Four_In_One!` evaluates the pressure-gradient
force, surface-pressure tendency, $M$, and $\omega$ in one column loop so that
they use the same discrete mass convergence.

### Thermodynamic equation

The adiabatic thermodynamic tendency is

```math
\frac{\partial T}{\partial t}
=-\mathbf{v}\cdot\nabla_h T
\mathcal{V}(T)
+\kappa T_v\frac{\omega}{p},
```

where $\mathcal{V}$ is the centered pressure-coordinate vertical-advection
operator and $\kappa=R_d/c_p$. Physical heating is not inserted into this
dynamical right-hand side: parameterizations update the provisional next state
after the dynamics step, as described in [Physics--dynamics coupling](@ref).

Hydrostatic geopotential is integrated upward from the surface geopotential:

```math
\Phi_{k-1/2}=\Phi_{k+1/2}
+R_dT_{v,k}\left(\ln p_{k+1/2}-\ln p_{k-1/2}\right),
```

and the full-level value is

```math
\Phi_k=\Phi_{k+1/2}
+R_dT_{v,k}\left(\ln p_{k+1/2}-\ln p_k\right).
```

Full-level log pressure is the Simmons--Burridge layer-mean value rather than
the logarithm of the arithmetic midpoint. At a zero-pressure model top the
top-layer limit is evaluated analytically.

### Moisture transport

Specific humidity is deliberately grid-only. It is advanced over the same
effective interval as the dry leapfrog state, starting from the current field
during the startup step and from the previous field during mature leapfrog
steps.

`Advance_Grid_Tracer!` applies two sequential finite-volume stages:

1. multidimensional spherical Van Leer transport in longitude and latitude;
2. monotone, nonuniform-grid piecewise parabolic transport through the vertical
   interfaces using $M$.

Both stages use a donor-cell positive baseline and limit only the
antidiffusive high-order flux. They subcycle automatically when the incoming
CFL bound requires it, allow a swept flux to cross complete cells, remove only
roundoff-scale undershoots, and fail loudly for a genuine negative result.
After transport, a global multiplicative restoration can compensate for
roundoff drift without changing the field's shape.

## Horizontal discretization

### Spherical harmonics and triangular storage

A scalar field $\chi(\lambda,\mu)$ is represented by

```math
\chi(\lambda,\mu)
=\sum_{m=0}^{M}\sum_{n=m}^{N}
 \widehat{\chi}_n^m P_n^m(\mu)e^{im\lambda},
\qquad \mu=\sin\phi.
```

Only nonnegative $m$ need to be stored because grid fields are real. The driver
constructs `num_spherical = num_fourier + 1`; the extra total-wavenumber row is
needed by vector transforms and derivative recurrences. Entries with $n<m$
are unused.

The horizontal Laplacian acts exactly on each retained scalar mode:

```math
\nabla_h^2\widehat{\chi}_n^m
=-\frac{n(n+1)}{a^2}\widehat{\chi}_n^m
=\Lambda_n\widehat{\chi}_n^m.
```

### Gaussian transform grid

The physical grid contains `nλ = 2nθ` uniformly spaced longitudes. The
latitudes are the roots of the degree-`nθ` Legendre polynomial in
$\mu=\sin\phi$, with Gaussian quadrature weights. The driver enforces

```math
n_\phi \ge \left\lceil\frac{3M+1}{2}\right\rceil,
```

which is the usual transform-grid requirement for quadratic nonlinear terms.

Synthesis uses an associated-Legendre transform followed by an inverse FFT;
analysis uses an FFT followed by Gaussian-weighted Legendre projection.
Nonlinear products, parameterizations, and diagnostics are evaluated on this
grid. Derivatives, Laplacians, and the inversion between $(\zeta,\delta)$ and
$(u,v)$ use precomputed spectral operators.

## Vertical discretization

The primitive-equation mode uses a Lorenz grid: $u$, $v$, $T$, $q$, $p$, and
$\omega$ are on full levels; $p_{k+1/2}$ and $M_{k+1/2}$ are on interfaces.
Level index 1 is at the model top and level `nd` is next to the surface.

### Coordinate choices

`vert_coord_option` selects the interface coefficients:

| Value | Description | Restrictions |
|:---|:---|:---|
| `"even_sigma"` | $A=0$ and evenly spaced $B\in[0,1]$ | arbitrary `nd` |
| `"uneven_sigma"` | stretched pure-sigma grid controlled internally by `scale_heights`, `surf_res`, and `exponent` | arbitrary `nd` |
| `"hybrid"` | smooth sine-squared blend from fixed pressure aloft to sigma near the surface | arbitrary `nd` |
| `"simmons_and_burridge"` | fixed 20-layer coefficient table | requires `nd = 20` |
| `"mcm"` | fixed legacy model coefficient table | requires the table's native level count |
| `"v197"` | fixed legacy V197 sigma table | requires the table's native level count |

The current driver exposes only the option name through `Model_Config`; the
optional stretching arguments of the lower-level `Vert_Coordinate` constructor
retain their source defaults.

`vert_difference_option` must currently be `"simmons_and_burridge"` for the
primitive-equation pressure and geopotential kernels. Despite its historical
configuration name, `vert_ref_level_option` is passed to the vertical advection
operator. Supported values are:

- `"second_centered_wts"`: second-order interface interpolation weighted by
  the two neighboring $\Delta p$ values;
- `"second_centered"`: arithmetic averaging of the neighboring full levels.

Both use no-flux top and bottom boundaries.

## Time integration

### Startup and leapfrog steps

The first call uses a one-step startup of length $\Delta t$: explicit terms use
the current state and the implicit coupling uses $\xi=\alpha\Delta t$.
Subsequent calls use the three-level leapfrog update over the interval
$2\Delta t$:

```math
X^{n+1}=X^{n-1}+2\Delta t\,\dot X^n.
```

The current level is filtered with coefficient $\nu_{RA}$:

```math
X^n\leftarrow X^n
+\nu_{RA}\left(X^{n-1}-2X^n+X^{n+1}\right).
```

The implementation calls this a filtered leapfrog scheme. The code does not
apply the Williams correction, so `robert_coef` should be interpreted as a
Robert--Asselin coefficient.

### Semi-implicit gravity-wave solve

Only the primitive-equation divergence--temperature--surface-pressure coupling
is treated semi-implicitly. The solver linearizes about
$p_{s,r}=10^5$ Pa and the constant reference profile $T_r=300$ K. The
hydrostatic, thermodynamic, and continuity linearizations form a vertical
matrix $\mathbf{B}$. For every total wavenumber $n$, the corrected divergence
tendency is obtained from

```math
\left[\mathbf{I}-\xi^2\Lambda_n\mathbf{B}\right]
\widehat{\dot{\boldsymbol\delta}}_n^m
=\widehat{\mathbf{R}}_n^m,
```

where $\xi=\alpha\Delta t$ during startup and
$\xi=2\alpha\Delta t$ afterward. Because $\Lambda_n$ depends only on $n$, one
small vertical inverse is precomputed per total wavenumber and shared by all
zonal wavenumbers $m$. Temperature and log-surface-pressure tendencies are
then recovered by back-substitution.

The shallow-water mode uses the analogous scalar Helmholtz correction for the
linear $-\nabla_h^2 h$ and $-h_0\delta$ gravity-wave terms. The barotropic mode
has no semi-implicit solve.

### Spectral damping

Vorticity, divergence, and temperature (or shallow-water height) receive an
implicit scale-selective damping. If $r$ is `damping_order` and $K_d$ is
`damping_coef`, the stored modal rate is

```math
d_n=K_d\left(\frac{\Lambda_n}{\Lambda_{N}}\right)^r.
```

The constructor requires even $r$. Damping is normalized by the largest
retained Laplacian eigenvalue, so `damping_coef` acts as the grid-scale rate.
The damping is included in the tendency through an implicit denominator,
which avoids a separate explicit diffusion stability limit.

## Global correction operators

The primitive-equation path can apply three global corrections after its
dynamics update:

- surface pressure is multiplied by a constant to restore global mean $p_s$;
- a spatially uniform temperature increment restores the pressure-mass-weighted
  integral of $c_pT+K$;
- nonnegative grid humidity is multiplied by a constant to restore its
  pressure-mass-weighted integral.

The targets come from the unfiltered previous physical state. These operators
compensate for spectral projection, filtering, and roundoff; they are numerical
corrections, not physical parameterizations. They are controlled by
`do_mass_correction`, `do_energy_correction`, and `do_water_correction`.

## Barotropic model

The barotropic mode advances only relative vorticity. With nondivergent wind
diagnosed from $\zeta$, its resolved equation is

```math
\frac{\partial\zeta}{\partial t}
=-\mathbf{v}\cdot\nabla_h(\zeta+f).
```

The code evaluates the vector-invariant flux $(\eta v,-\eta u)$ on the
Gaussian grid, takes its rotational part, applies spectral damping, advances
with filtered leapfrog, and reconstructs a zero-divergence wind.

## Shallow-water model

The one-layer model advances $\zeta$, $\delta$, and geopotential height $h$:

```math
\frac{\partial\mathbf{v}}{\partial t}
=-\eta\mathbf{k}\times\mathbf{v}-\nabla_h(K+h)
-\kappa_m\mathbf{v},
```

```math
\frac{\partial h}{\partial t}
=-\mathbf{v}\cdot\nabla_h h-h\delta
-\kappa_t(h-h_{eq}).
```

Here $h$ is stored internally in fields whose primitive-equation names contain
`ps` or `lnps`; the output mapping exposes it as `:h`. In the numerical
equations $h$ is added to $K$, so it behaves as layer geopotential with units
m$^2$ s$^{-2}$. The current NetCDF metadata instead labels it geopotential
height in metres; users should account for this metadata inconsistency when
analyzing shallow-water output. The linear drag and height relaxation are
intrinsic forcing for this idealized equation set, not part of the
primitive-equation physics suite. Their parameters are `kappa_m`, `kappa_t`,
and `h_0`.

## Primitive-equation step sequence

One call to `Atmosphere_Update!` performs the following operations:

1. diagnose pressure, $T_v$, geopotential, $\omega$, and interface mass flux;
2. assemble horizontal and vertical adiabatic tendencies on the grid;
3. transport grid humidity over the effective leapfrog interval;
4. transform dry tendencies to spectral space and construct vorticity and
   divergence tendencies;
5. solve the semi-implicit Helmholtz problem;
6. apply spectral damping and filtered leapfrog updates;
7. apply the enabled dynamics conservation corrections;
8. apply the ordered gridpoint parameterizations to the provisional next
   state;
9. project the post-physics dry fields back to spectral space, restore any
   projection drift in the post-physics mass and energy integrals, preserve the
   post-physics moisture integral, and rotate time levels.

This sequence is the boundary between the dynamical core and the
parameterizations. It also explains why physics tendencies are not added to a
leapfrog right-hand side.
