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
temperature and humidity relaxation rates. The reference is Frierson's
Simple Betts--Miller formulation (2007), with the O'Gorman and Schneider (2008)
modifications implemented in Isca's `qe_moist_convection.F90`: consistent
virtual-temperature buoyancy and exact vapor-pressure/specific-humidity
conversion. It uses Isca's shallower shallow-convection option.

### Parcel ascent and triggering

For each column, pressure increases from model top to surface. The surface
parcel starts with the lowest-level $T$ and mixing ratio $r=q/(1-q)$.

- A supersaturated starting parcel is adjusted to saturation at the surface.
- An unsaturated parcel follows a dry adiabat to its lifting condensation
  level (LCL), located by bisection in log pressure. Below the LCL its actual
  mixing ratio remains equal to the starting value; the saturation mixing
  ratio used for the reference profile is stored separately. The LCL may be
  at the first full level. If the parcel remains unsaturated through the top
  full level, it retains its dry-adiabatic profile and receives no moist
  adjustment (`lcl = 0`).
- Above the LCL, a saturated moist adiabat is integrated with a second-order
  Runge--Kutta step in log pressure. Saturation at the RK midpoint is evaluated
  at the arithmetic midpoint pressure.
- With `use_virtual_temperature = true`, parcel and environmental buoyancy
  use $T_v=T[1+(R_v/R_d-1)q]$, consistent with the dynamics. With the flag
  disabled, buoyancy uses actual temperature. The ascent and enthalpy
  calculation always use actual temperature.
- Discrete buoyancy is integrated using $R_d(T_{v,p}-T_v)\,\Delta\ln p$ to diagnose
  convective inhibition and CAPE. The first buoyant level is the level of free
  convection (LFC); the first stable level above a contiguous buoyant region
  terminates it. If buoyancy reaches the model top, the top full level is the
  level of zero buoyancy (LZB). Roundoff-scale differences are treated as
  neutral, and neutral levels do not establish an LFC.

A column is inactive if it has no positive contiguous CAPE, if the initially
dry parcel has no water vapor, or if the parcel becomes colder than 173.16 K
before reaching buoyancy. `Betts_Miller_Column` returns actual
`parcel_mixing_ratio` and separate `parcel_saturation_mixing_ratio`
diagnostics. Above the end of ascent, temperature and actual mixing ratio
retain their input values; saturation mixing ratio is zero at unvisited levels.

The returned `cin` follows Isca's signed integration below an in-domain LCL;
buoyant dry layers can reduce it. No-CAPE and pre-buoyancy cold-cutoff paths
reset CIN to zero, while a parcel with no represented LCL retains the positive-only dry
inhibition diagnostic. CIN is not an adjustment threshold and should not be
interpreted as a universal positive-only column integral.

### Reference state and relaxation

From the LZB through the surface, the parcel temperature is the preliminary
reference temperature. The preliminary reference humidity uses exact relative
humidity $\mathcal{H}_{BM}$ by scaling saturation vapor pressure:

```math
T_{ref,k}=T_{p,k},\qquad
e_{ref,k}=\mathcal{H}_{BM}e_s(T_{p,k}),\qquad
q_{ref,k}=\frac{\epsilon e_{ref,k}}
{p_k-(1-\epsilon)e_{ref,k}},\qquad \epsilon=\frac{R_d}{R_v}.
```

This definition also applies below the LCL: reference humidity depends on
saturation at the dry parcel temperature, while the actual parcel mixing
ratio remains conserved. Multiplying saturation mixing ratio by
$\mathcal{H}_{BM}$ would only approximate the requested relative humidity.
Above the adjustment layer, reference temperature and humidity retain their
environmental values.

The unbalanced relaxation rates are

```math
\left.\frac{\partial T_k}{\partial t}\right|_{BM}
=\frac{T_{ref,k}-T_k}{\tau_{BM}},
\qquad
\left.\frac{\partial q_k}{\partial t}\right|_{BM}
=\frac{q_{ref,k}-q_k}{\tau_{BM}}.
```

### Deep-convection enthalpy closure

With layer mass $m_k=(p_{k+1/2}-p_{k-1/2})/g$, the preliminary rates give
$P_q=-\sum_k\dot q_k m_k$ and $P_T=(c_p/L_v)\sum_k\dot T_k m_k$.
Deep convection requires both integrals to be positive. Roundoff tolerances
scale with the sum of the absolute terms; no physical CAPE cutoff is imposed.

The default `bm_energy_correction = :isca` follows Isca's Simple Betts--Miller
deep closure. If $P_q>P_T$, humidity rates are multiplied by $P_T/P_q$ and
precipitation is $P_T$. Otherwise, precipitation is $P_q$ and a uniform
temperature-rate correction is applied over the adjustment layer:

```math
C=-\frac{\sum_k(c_p\dot T_k+L_v\dot q_k)m_k}{c_p\sum_k m_k}.
```

The temperature reference receives the corresponding $\tau_{BM}C$ shift.
Reference humidity continues to use the preliminary parcel temperature; it is
not recomputed after the shift. Selecting `:timescale` retains the previous
deep closure: scale whichever precipitation-equivalent adjustment is larger
down to the smaller. Both choices conserve moist enthalpy and report moisture
loss as precipitation, but can produce different heating profiles at the same
rainfall. Scaled rates correspond to a longer effective relaxation time toward
the original reference, rather than a new reference at the nominal time scale.

### Conservative shallow convection

Positive CAPE and $P_T>0$ with $P_q\le0$ can produce zero-rain shallow
transport. Starting at the diagnosed LZB, upper layers are excluded until
moistening in a transition cell can balance drying in the fully included
lower cells. The included fraction is

```math
f Q_{top}+Q_{lower}=0,\qquad
Q_{top}=\dot q_{top}m_{top},\qquad 0\le f\le1.
```

Both preliminary cell rates are multiplied by $f$ in the transition cell;
rates above it are zero. A uniform temperature correction then conserves
moist enthalpy over the remaining full grid cells. The fraction enters the
rates once, and is not applied again to the correction mass. Consequently,
shallow transport has zero integrated moisture and thermal tendency and
exactly zero reported rain. Both deep energy-correction options use this
same shallow closure.

An exact or roundoff-scale zero $P_q$ keeps the full feasible adjustment
region. Exact interface crossings exclude the zero-fraction cell. If there
is no feasible depth, or only one cell remains, the column receives no
adjustment. Nonpositive $P_T$ also produces no adjustment. The Isca closure
approaches the zero-drying boundary continuously. The alternative timescale
closure can change its heating profile abruptly there because its deep
thermal rates approach zero, while shallow convection retains conservative
heat redistribution.

For shallow convection, `reference_temperature` and `reference_humidity`
are effective cell targets satisfying $X_{ref}=X+\tau_{BM}\dot X$ after
fractional penetration and correction. They can differ from the preliminary
parcel/RH targets. Inactive columns retain environmental references.

`Betts_Miller_Column` reports `regime` as `:none`, `:deep`, or `:shallow`.
`active` includes nonzero zero-rain shallow adjustment. `lzb` remains the
diagnosed buoyancy limit; `adjustment_top` is the first adjusted full level
and `top_fraction` is its included preliminary fraction. For deep adjustment,
these are the LZB and 1; for inactive columns, they are 0 and 0. The corrected
temperature rate in a fractional cell includes the uniform correction as
well as its fraction of the preliminary rate.

Physics applies these rates explicitly for one $\Delta t_p$ substep. The driver
therefore requires `config.Δt <= bm_tau`, bounding the explicit relaxation
fraction by one. This bound does not by itself establish positive temperature
or stability on arbitrary vertical grids. The ascent follows Isca's RK2
discretization; unusually coarse pressure grids require an accuracy check.

| Key | Default | Meaning |
|:---|:---|:---|
| `"do_Betts_Miller"` | `false` | enable convective adjustment |
| `"bm_tau"` | `7200.0` | relaxation time $\tau_{BM}$, s |
| `"bm_relative_humidity"` | `0.8` | reference relative humidity $\mathcal{H}_{BM}\in(0,1]$ |
| `"bm_energy_correction"` | `:isca` | deep enthalpy closure: `:isca` or `:timescale`; strings accepted |

The optional `"initial_humidity_floor"` belongs to the `:Moist_Spinup`
initial condition, not to the convection calculation. It can suppress
very dry represented points. The current initial humidity is a grid tracer
and is not spectrally projected.

Physical column inputs have positive temperature, increasing full/interface
pressures with each full level inside its interfaces, and `0 <= q < 1`.
Standalone column calls retain legacy handling of finite negative humidity:
they calculate relative to `max(q,0)` without mutating the input. Conservation
identities then refer to that cleaned humidity. Material negative vapor mass
is not a physically valid input; coupled physics validates its updated state.

### Validation and migration

The implementation conserves pressure-mass-weighted $c_pT+L_vq$ during BM
adjustment and accounts for water loss as precipitation. These are fixed-mass
column identities. A model budget also includes surface exchange, Newtonian
relaxation, pressure adjustment, PBL mixing, dynamics, and spectral/numerical
corrections; it cannot be inferred from the BM identity alone.

Independent fixtures and an optional Fortran harness cover 87 Isca columns,
including both deep branches and shallow transport. Comparing with a control
that changes only Isca's LCL lookup to bisection gives temperature-rate
agreement within `2e-17 K/s`. Differences from the original lookup and its
boundary behavior are recorded explicitly. See
[Betts--Miller validation](betts-miller-validation.md) for the static audit,
five-day T21 timestep comparison, budget scope, and reproduction commands.

Corrected virtual buoyancy, exact reference RH, and shallow transport change
the model's results. Existing experiments need rerunning. Selecting
`:timescale` retains the previous deep closure, but does not recover the old
scheme's complete behavior. Short integration checks establish stability for
their prescribed case; they do not establish a climatology or tune RH/tau.

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

The optional LRF adds a prescribed longwave temperature response to humidity
anomalies. It changes temperature only and uses fixed coefficients and
reference humidity throughout a run. Three artifact schemes are supported.
Each is evaluated independently in every model column, with matrix rows
representing temperature response levels and columns representing humidity
perturbation levels. There is no horizontal convolution or online zonal
averaging.

### Legacy linear humidity response

The `linear_q_v1` scheme uses `LRF_State`. For each latitude and longitude,

```math
\left.\frac{\partial T_{k_o}}{\partial t}\right|_{LRF}
=\frac{1}{t_{day}}
\sum_{k_i=1}^{N}
L_{k_o k_i}(\phi)
\left[q_{k_i}-q_{ref,k_i}(\lambda,\phi)\right].
```

The matrix values have units K day$^{-1}$ per unit specific humidity in kg/kg.
A legacy JLD2 file contains:

- `LRF_LW_q` with size `(nd, nd, nθ)`;
- `ref_q` with size `(nλ, nθ, nd)`.

A file with no `scheme` key is interpreted as `linear_q_v1`, preserving
compatibility with existing artifacts.

### Regularized logarithmic response

Both modified schemes define a natural-log humidity anomaly

```math
\delta\chi_k(\lambda,\phi)
=\ln\left[q_k(\lambda,\phi)+q_0\right]-\chi_{ref,k}(\phi).
```

Here $q$ and $q_0$ are numerical values of specific humidity in kg/kg. The
positive $q_0$ makes the logarithm finite at zero humidity. It does not alter
the prognostic humidity or limit heating. Negative and nonfinite humidity
inputs are rejected. The response is linear in $\delta\chi$, while its local
derivative with respect to $q$ scales as $1/(q+q_0)$.

The shared-matrix scheme, `regularized_log_tapered_v1`, uses
`Regularized_LRF_State` and computes

```math
\left.\frac{\partial T_{k_o}}{\partial t}\right|_{LRF}
=\frac{\alpha\,w(\phi)}{t_{day}}
\sum_{k_i=1}^{N} B^*_{k_o k_i}\,\delta\chi_{k_i}(\lambda,\phi).
```

The stored `taper` supplies $w(\phi)$ in $[0,1]$; a zero weight gives exactly
zero heating at that latitude. The runtime does not generate a taper from
latitude or assume a particular tropical transition.

The latitude-specific scheme, `regularized_log_latitude_v1`, uses
`Latitude_LRF_State` and computes

```math
\left.\frac{\partial T_{k_o}}{\partial t}\right|_{LRF}
=\frac{\alpha}{t_{day}}
\sum_{k_i=1}^{N} B_{k_o k_i}(\phi)\,\delta\chi_{k_i}(\lambda,\phi).
```

Every latitude has a separate vertical matrix, shared by all its longitudes.
This scheme has no separate taper. Both logarithmic matrices have units
K day$^{-1}$ per unit log-humidity anomaly; all schemes divide by
`config.day_to_sec` to return K s$^{-1}$.

The logarithmic JLD2 artifact layouts are:

| Key | `regularized_log_tapered_v1` | `regularized_log_latitude_v1` |
|:---|:---|:---|
| `scheme` | `"regularized_log_tapered_v1"` | `"regularized_log_latitude_v1"` |
| Response matrix | `B_star(nd, nd)` | `B_by_lat(nd, nd, nθ)` |
| `chi_reference` | `(nd, nθ)` | `(nd, nθ)` |
| `taper` | vector of length `nθ`, in `[0, 1]` | unused |
| `q0_kg_kg` | positive finite scalar | positive finite scalar |
| `alpha` | nonnegative finite scalar | nonnegative finite scalar |
| `latitude` | vector of length `nθ`, in degrees | vector of length `nθ`, in degrees |

Matrices and references must be finite. `alpha` and `q0_kg_kg` come from the
artifact; the driver has no strength override or activation ramp. Model
indices, including vertical ordering, must match the artifact: the loader
does no regridding. On the driver path, stored latitude values and ordering
must match `rad2deg.(mesh.θc)` within `1e-10` degrees. Direct calls to
`Load_LRF_State(...; latitude=nothing)` skip this coordinate check. Legacy
artifacts are dimension-checked but do not undergo the logarithmic schemes'
finite-value or latitude-coordinate checks. Stored pressure coordinates and
provenance metadata are not validated by the runtime loader.

For a reference centered on a calibration climate, use
$\chi_{ref}=\langle\ln(q+q_0)\rangle$ over the chosen time and longitude
sample. This generally differs from $\ln(\langle q\rangle+q_0)$.
With fixed matrices, the former makes mean added heating vanish over that
sample at each level and latitude. A different climate or reference policy
can have nonzero mean heating. If $L$ is an offline humidity Jacobian, the
transformed matrix can be formed as
$B_{k_o k_i}=L_{k_o k_i}(q_{cal,k_i}+q_0)$; the model consumes the resulting
matrix and does not run a radiation model online.

### Configuration and coupling

Enable LRF in `Model_Config` with `moisture_processes = true` and

```julia
physics_params = Dict{String,Any}(
    "do_LRF" => true,
    "LRF_file" => "/path/to/latitude_lrf.jld2",
)
```

The driver loads and validates the file once before integration. A minimal
latitude-specific artifact can be written from matching model-grid arrays:

```julia
using JLD2

JLD2.jldsave("latitude_lrf.jld2";
    scheme = "regularized_log_latitude_v1",
    B_by_lat = B_by_lat,
    chi_reference = chi_reference,
    q0_kg_kg = 1.0e-8,
    alpha = 1.0,
    latitude = rad2deg.(mesh.θc),
)
```

LRF sees humidity after convection, condensation, surface exchange, and PBL
mixing in each physics substep. After Held--Suarez forcing it adds the
explicit increment `T += physics_dt * lrf_tendency`, before physical-state
validation and dry-air adjustment. The later spectral synchronization uses
the post-physics energy target, retaining this physical heating. The output
symbol `:lrf_dt` writes `lrf_dta_dt` in K s$^{-1}$: the physics-interval mean
of the substep rates, subsequently averaged over the output interval. Its
buffer is zero when LRF is disabled.

Coefficients are reloaded from the configured external artifact on a restart;
they are not stored in checkpoints. Keep the artifact and its provenance
with the experiment. These prescribed humidity responses do not include an
online cloud, shortwave, or temperature-radiation calculation, and finite
zero-humidity heating does not establish stability or accuracy far from the
calibration climate.

`post_processing/Build_Experiment_LRF.py` prepares separately centered
T42L20 artifacts for `ctrl_BM`, `sst1.0_BM`, and `sst2.5_BM` from native
archives. It requires a climlab/RRTMG environment and the external utilities
selected by `--utilities`; those utilities are not bundled with JGCM.
`post_processing/Write_Experiment_LRF.jl STAGING_DIRECTORY` exports its binary
staging inputs to JLD2 and checks the actual model grid, Python/Julia heating
agreement, and reference/zero-humidity behavior. Offline HDF5 calibration
files must be exported to JLD2 before using them as `LRF_file`.

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
