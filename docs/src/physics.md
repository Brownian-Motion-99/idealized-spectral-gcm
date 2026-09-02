# Physical parameterizations

This page describes only the parameterized processes implemented in
`src/Physics`: Held--Suarez forcing, Betts--Miller convection, large-scale
condensation, surface exchange, boundary-layer mixing, and the moisture
linear-response forcing. The resolved equations and numerical core are in
[Dynamical core](@ref). Symbols follow [Notation](@ref).

The parameterizations are available only to `model_type = :PrimitiveEquation`.
Their switches and coefficients live in `Model_Config.physics_params`.

## Physics--dynamics coupling

The dynamical core first constructs a provisional next state. Physics then
updates a private Gaussian-grid copy of $(u,v,T,q)$ directly. It does **not**
add parameterized rates to the leapfrog right-hand side.

The effective dynamics interval is $\Delta t$ during the startup step and
$2\Delta t$ during a mature leapfrog step. Physics always uses the base model
timestep `config.Δt`, denoted $\Delta t_p$, so a mature leapfrog interval
contains two ordered physics substeps. The process order in every substep is:

1. Betts--Miller convection;
2. large-scale condensation;
3. surface sensible heating;
4. surface evaporation;
5. implicit boundary-layer mixing;
6. Held--Suarez Rayleigh friction, including frictional heating;
7. Held--Suarez Newtonian relaxation;
8. moisture linear-response-function heating.

This order is observable. For example, saturation adjustment sees the
post-convection column, and Held--Suarez thermal relaxation sees frictional
heating and all moist-process temperature changes from the same substep.
Tendency and flux diagnostics are averaged over the substeps, rather than
reporting only the last substep.

After each moist substep, column surface pressure and humidity are adjusted so
that dry-air mass is unchanged and the water amount produced by physics is
retained in every layer. After all physics, dry fields are projected back to
the spectral truncation and humidity remains grid-only.

## Shared moist thermodynamics

Betts--Miller convection, large-scale condensation, and surface evaporation use
the same saturation functions. Let $e_s(T)$ be the Smithsonian saturation
vapor pressure, evaluated over ice below 253.16 K, over liquid above 273.16 K,
and with a linear phase blend between those temperatures. With
$\epsilon=R_d/R_v$, the exact saturation mixing ratio and saturation specific
humidity are

```math
r_s(T,p)=\frac{\epsilon e_s(T)}{p-e_s(T)},
```

```math
q_s(T,p)=\frac{\epsilon e_s(T)}
{p-(1-\epsilon)e_s(T)}.
```

$r_s$ is water-vapor mass per unit dry-air mass and is used in the parcel
calculation. $q_s$ is water-vapor mass per unit moist-air mass and is used for
the prognostic humidity. Keeping these two quantities distinct avoids the
small but systematic error caused by treating mixing ratio as specific
humidity.

## Held--Suarez forcing

Enable both the thermal and momentum parts with `"do_HS_Forcing" => true`.
Configured rates `k_a`, `k_s`, and `k_f` are in day$^{-1}$ and are divided by
`day_to_sec` internally.

### Newtonian thermal relaxation

The equilibrium temperature is

```math
T_{eq}(\phi,p)=\max\left\{
T_{strat},
\left[T_{eq,0}-\Delta T_y\sin^2\phi
-\Delta\theta_z\cos^2\phi\ln\left(\frac{p}{p_0}\right)\right]
\left(\frac{p}{p_0}\right)^\kappa
\right\},
```

where $p_0=10^5$ Pa. The relaxation rate is

```math
k_T(\phi,\sigma)=k_a+
(k_s-k_a)\max\left(0,\frac{\sigma-\sigma_b}{1-\sigma_b}\right)
\cos^4\phi.
```

The explicit finite update is

```math
T^{new}=T^{old}+\Delta t_p k_T(T_{eq}-T^{old}).
```

The code requires $0\leq\Delta t_p k_T\leq1$. `grid_t_eq`, output as `:t_eq`,
contains the diagnosed equilibrium temperature.

### Rayleigh friction and frictional heating

The low-level wind damping rate is

```math
k_v(\sigma)=k_f\max\left(0,
\frac{\sigma-\sigma_b}{1-\sigma_b}\right).
```

The wind update is explicit:

```math
(u^{new},v^{new})=(1-\Delta t_p k_v)(u^{old},v^{old}),
```

and requires $\Delta t_p k_v\leq1$. The exact kinetic-energy loss of this
finite update is returned locally to temperature:

```math
T^{new}\leftarrow T^{new}
+\frac{(u^{old})^2+(v^{old})^2-(u^{new})^2-(v^{new})^2}{2c_p}.
```

Thus the friction step conserves local kinetic plus sensible energy to
roundoff; the subsequent Newtonian relaxation can still add or remove energy.

### Held--Suarez parameters

| Key | Example value | Meaning |
|:---|:---|:---|
| `"do_HS_Forcing"` | `true` | enable Rayleigh friction and Newtonian relaxation |
| `"σ_b"` | `0.7` | sigma coordinate at the top of the frictional layer |
| `"k_f"` | `1.0` | surface momentum damping rate, day$^{-1}$ |
| `"k_a"` | `1/40` | free-atmosphere thermal damping rate, day$^{-1}$ |
| `"k_s"` | `1/4` | near-surface thermal damping rate, day$^{-1}$ |
| `"T_equator"` | `294.0` | $T_{eq,0}$, K |
| `"T_stratosphere"` | `200.0` | lower bound $T_{strat}$, K |
| `"ΔT_y"` | `60.0` or `65.0` | equator-to-pole contrast, K |
| `"Δθ_z"` | `10.0` | vertical potential-temperature contrast, K |

The six coefficients without an internal `get(..., default)` call must be
present whenever the scheme is enabled. Use the Unicode keys shown above;
ASCII aliases such as `"sigma_b"` are not read by the forcing routine.

## Betts--Miller convection

The Betts--Miller implementation diagnoses a lifted surface parcel, identifies
a contiguous buoyant layer, constructs reference profiles, and returns
temperature and humidity relaxation rates.

### Parcel ascent and triggering

For each column, pressure increases from model top to surface. The surface
parcel starts with the lowest-level $T$ and mixing ratio $r=q/(1-q)$.

- A supersaturated starting parcel is adjusted to saturation at the surface.
- An unsaturated parcel follows a dry adiabat to its lifting condensation
  level (LCL), located by bisection in log pressure.
- Above the LCL, a saturated moist adiabat is integrated with a second-order
  Runge--Kutta step in log pressure. Saturation at the RK midpoint is evaluated
  at the arithmetic midpoint pressure.
- Discrete buoyancy is integrated using $R_d(T_p-T)\,\Delta\ln p$ to diagnose
  convective inhibition and CAPE. The first buoyant level is the level of free
  convection (LFC); the first stable level above a contiguous buoyant region
  terminates it. If buoyancy reaches the model top, the top full level is the
  level of zero buoyancy (LZB).

A column is inactive if it has no positive contiguous CAPE, if the initially
dry parcel has no water vapor, or if the parcel becomes colder than 173.16 K
before reaching buoyancy.

### Reference state and relaxation

From the LZB through the surface, the parcel temperature is the reference
temperature and the reference mixing ratio is a fixed fraction
$\mathcal{H}_{BM}$ of parcel saturation:

```math
T_{ref,k}=T_{p,k},\qquad
r_{ref,k}=\mathcal{H}_{BM}r_{s,k},\qquad
q_{ref,k}=\frac{r_{ref,k}}{1+r_{ref,k}}.
```

The unbalanced relaxation rates are

```math
\left.\frac{\partial T_k}{\partial t}\right|_{BM}
=\frac{T_{ref,k}-T_k}{\tau_{BM}},
\qquad
\left.\frac{\partial q_k}{\partial t}\right|_{BM}
=\frac{q_{ref,k}-q_k}{\tau_{BM}}.
```

The column-integrated precipitation-equivalent rates inferred from drying and
heating are compared. The larger adjustment is scaled down so that the two
energy-equivalent precipitation rates agree. A column is rejected if either
unbalanced integral is nonpositive. The accepted common value is reported as
the convective precipitation flux.

Physics applies these rates explicitly for one $\Delta t_p$ substep. The driver
therefore requires `config.Δt <= bm_tau`, preventing an individual relaxation
step from passing its reference profile.

| Key | Default | Meaning |
|:---|:---|:---|
| `"do_Betts_Miller"` | `false` | enable convective adjustment |
| `"bm_tau"` | `7200.0` | relaxation time $\tau_{BM}$, s |
| `"bm_relative_humidity"` | `0.8` | reference relative humidity $\mathcal{H}_{BM}\in(0,1]$ |

The optional `"initial_humidity_floor"` belongs to the `:Moist_Spinup`
initial condition, not to the convection calculation. It can suppress
roundoff-scale dry points in a spectrally truncated analytic initialization.

## Large-scale condensation

Enable saturation adjustment with `"do_Lscale_Cond" => true`. At every
supersaturated grid cell, the scheme linearizes $q_s$ about the current
$(T^*,p)$ and computes

```math
\Delta q_{LS}=
\frac{q_s(T^*,p)-q^*}
{1+\mathcal{L}(L_v/c_p)\left.\partial q_s/\partial T\right|_{T^*,p}},
```

```math
\Delta T_{LS}=-\mathcal{L}\frac{L_v}{c_p}\Delta q_{LS}.
```

The heating fraction $\mathcal{L}$ is either a scalar or an `nλ × nθ` array,
and every value must lie in $[0,1]$. It is read first from
`"condensation_heating_fraction"`; the legacy key `"L"` is used as a fallback,
with a final default of 1.

The scheme updates $T$ and $q$ immediately. Removed vapor becomes positive
precipitation,

```math
P_{LS}=-\frac{1}{\Delta t_p}
\sum_k \Delta q_{LS,k}\frac{\Delta p_k}{g},
```

and the diagnostic liquid-water-content field is the positive condensation
rate $-\Delta q_{LS}/\Delta t_p$. Subsaturated cells are unchanged. When
$\mathcal{L}<1$, both the temperature increment and the saturation-adjustment
denominator use the reduced heating, so exact moist-static-energy closure is
not intended.

When Betts--Miller is also enabled, convection is applied first and large-scale
condensation sees its updated state. Total precipitation is the sum of the two
schemes; `:bm_precip` isolates the convective contribution.

## Surface exchange

Surface sensible heat and evaporation use the lowest full level, the local
wind speed $V_c=(u^2+v^2)^{1/2}$, and a hydrostatic lowest-level height proxy

```math
z_a=\frac{R_dT_v}{2g}
\ln\left(\frac{p_s}{p_{N-1/2}}\right).
```

The prescribed surface temperature is a callback
`lower_boundary_temperature(longitude, latitude)` with both coordinates in
radians. The default is zonally symmetric:

```math
T_s(\phi)=271+29\exp\left[-\frac{\phi^2}{2(26^\circ)^2}\right]\ \mathrm{K}.
```

### Sensible heat

With $\beta_H=C_HV_c\Delta t_p/z_a$, the backward-implicit bulk update is

```math
T_N^{new}=\frac{T_N^{old}+\beta_HT_s}{1+\beta_H}.
```

The reported upward surface sensible heat flux is diagnosed from the actual
finite layer-energy increment:

```math
SH=\frac{\Delta p_N}{g}\,c_p
\frac{T_N^{new}-T_N^{old}}{\Delta t_p}.
```

Enable it with `"do_Sensible_Heating" => true` and provide `"C_H"` (the
examples use `0.0044`).

### Evaporation

The saturated surface humidity is $q_s(T_s,p_s)$. A no-dew target ensures that
the surface can add water but cannot remove it:

```math
q_{target}=\max(q_N^{old},q_s(T_s,p_s)).
```

With $\beta_E=C_EV_c\Delta t_p/z_a$,

```math
q_N^{new}=\frac{q_N^{old}+\beta_Eq_{target}}{1+\beta_E}.
```

The upward latent heat flux is

```math
LH=\frac{\Delta p_N}{g}\,L_v
\frac{q_N^{new}-q_N^{old}}{\Delta t_p}.
```

Enable it with `"do_Surface_Evaporation" => true` and provide `"C_E"` (the
examples use `0.0044`). Evaporation is skipped when
`moisture_processes = false` even if its switch is true.

## Implicit boundary-layer mixing

The boundary-layer scheme vertically diffuses potential temperature and
specific humidity. It currently does **not** diffuse or drag momentum; the
Held--Suarez Rayleigh term is the available primitive-equation momentum drag.

The interface diffusivity is

```math
K_E=C_DV_cz_a
```

inside the prescribed boundary layer. Above it, $K_E$ decays as a Gaussian in
pressure with a hard-coded 10,000 Pa scale. The top can be selected in two
ways:

| Key combination | Interpretation |
|:---|:---|
| `"PBL_Top_Mode" => :PressureLevel`, `"PBL_Top_Value" => 85000.0` | constant mixing below a pressure interface in Pa |
| `"PBL_Top_Mode" => :ModelLevel`, `"PBL_Top_Value" => 4` | constant mixing through the lowest four model layers |

`PBL_Top_Value` must be a `Float64` in pressure mode and an `Int64` in model
level mode. Enable the scheme with `"do_Implicit_PBL_Scheme" => true` and
provide `"C_D"`.

For a transported scalar $\chi\in\{\theta,q\}$, the pressure-coordinate flux
coefficient is $g^2\rho^2K_E$. The backward-Euler finite-volume equation is

```math
\frac{\chi_k^{new}-\chi_k^{old}}{\Delta t_p}
=\frac{F_{k+1/2}^{new}-F_{k-1/2}^{new}}{\Delta p_k},
```

```math
F_{k+1/2}^{new}
=\frac{g^2\rho_{k+1/2}^2K_{E,k+1/2}}
{p_{k+1}-p_k}
\left(\chi_{k+1}^{new}-\chi_k^{new}\right).
```

Interface density uses pressure and the mean virtual temperature of the two
neighboring levels. The resulting tridiagonal column system is solved with a
Thomas algorithm. Temperature is converted to potential temperature for the
solve and converted back afterward.

## Moisture linear response function

The optional LRF represents a prescribed longwave temperature response to
humidity anomalies. For each latitude and longitude,

```math
\left.\frac{\partial T_{k_o}}{\partial t}\right|_{LRF}
=\frac{1}{t_{day}}
\sum_{k_i=1}^{N}
L_{k_o k_i}(\phi)
\left[q_{k_i}-q_{ref,k_i}(\lambda,\phi)\right].
```

The response matrix values are interpreted as K day$^{-1}$ per unit specific
humidity and converted using `day_to_sec`. The scheme changes temperature only.

Enable it with `"do_LRF" => true` and set `"LRF_file"` to a JLD2 file that
contains:

- `LRF_LW_q` with size `(nd, nd, nθ)`;
- `ref_q` with size `(nλ, nθ, nd)`.

The file is loaded and dimension-checked once before integration. LRF requires
`moisture_processes = true`; the resulting temperature rate is output with
`:lrf_dt`.

## Parameterization diagnostics

The parameterization-related `vars_to_output` symbols are:

| Symbol | Quantity | Units |
|:---|:---|:---|
| `:t_eq` | Held--Suarez equilibrium temperature | K |
| `:shflx` | upward surface sensible heat flux | W m$^{-2}$ |
| `:lhflx` | upward surface latent heat flux | W m$^{-2}$ |
| `:precip` | total precipitation flux | kg m$^{-2}$ s$^{-1}$ |
| `:bm_dt` | Betts--Miller temperature tendency | K s$^{-1}$ |
| `:bm_dq` | Betts--Miller humidity tendency | s$^{-1}$ |
| `:bm_precip` | Betts--Miller precipitation flux | kg m$^{-2}$ s$^{-1}$ |
| `:lrf_dt` | LRF temperature tendency | K s$^{-1}$ |

Disabled-process diagnostic buffers are reset to zero each dynamics call.
