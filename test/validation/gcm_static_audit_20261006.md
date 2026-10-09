# Whole-model static audit, 2026-10-06

This report describes the pre-repair October 6 source snapshot. F1 is now
corrected; see [the repair validation](f1_repair_20261006.md). The obsolete
baseline evidence script asserted faulty behavior and has been removed in
favor of the current post-physics regression tests.

The current GCM has confirmed state-coupling and diagnostic defects despite
passing its existing tests. The highest-priority issue is that post-physics
winds and spectral vorticity/divergence are synchronized, but grid vorticity
and divergence are not. Those stale grid fields enter the next dynamics step.
This affects the current Held–Suarez configuration, including `ctrl_BM`.

No additional error was established in the ordinary Betts–Miller deep/shallow
column closure. That narrower conclusion does **not** certify the complete
dynamics–physics system, and the newly identified coupling defect qualifies
the earlier interpretation of the Wheeler–Kiladis results. Static inspection
cannot establish how much of the ±5 signal it explains.

## Scope and evidence

Reviewed the current working tree based on commit `8ef276c`, including the
pre-existing changes to the LRF implementation, driver, experiment, and offline
interpolator. The review traced initialization, spherical transforms, primitive,
barotropic and shallow-water dynamics, pressure/hydrostatic calculations,
semi-implicit integration, tracer transport, global corrections, moist physics,
surface exchange/PBL, Held–Suarez forcing, LRF, restart, and NetCDF output.

Production source files and the existing user changes were left untouched.
The original audit added this report, a baseline evidence script, its results,
and a source-hash manifest. The experiments' long integrations were not launched.

Validation of that snapshot used Julia 1.8.5, the complete regression suite,
and the separate baseline evidence script.

The existing suite passes **15,626 assertions**: 90 core numerics, 170 runtime,
15,045 physics, and 321 integration assertions. The separate evidence script
passes **16 assertions that demonstrate the reported faulty behavior**; it is
not a regression suite asserting that the model is correct. It uses disposable
temporary files, small grids, analytic columns, and one-step integration.

See the [source manifest](gcm_static_audit_20261006_sources.toml). The observed
results are summarized here; the old failure-demonstration script and duplicate
console output are omitted from the repository.

Severity describes the consequence when the affected path is used. “Current
T42” refers to `exp/HSt42/HS.jl`: even sigma, 20 levels, 600 s timestep,
pressure-based PBL top at 850 hPa, 12-hour output, 50-day chunks, and built-in
moist initialization. LRF is conditional on the environment configuration.

| ID | Severity | Finding | Current T42 relevance |
|---|---|---|---|
| F1 | High | Grid vorticity/divergence remain stale after wind physics | Active with Held–Suarez drag |
| F2 | Medium | Vorticity/divergence tendency output is identically zero | If those output variables are requested |
| F3 | High for data integrity | Reusing an output directory can mix old and new runs | Conditional on directory reuse |
| F4 | Medium | Restart checks do not establish state/configuration compatibility | Warm starts, particularly changed timestep/grid |
| F5 | Medium | Initial-condition aliases produce inconsistent prognostics | Array input and built-in 2-D cases; not `:Moist_Spinup` |
| F6 | Medium | Energy fixer assumes the atmospheric top is at zero pressure | Finite-top grids; current even sigma unaffected |
| F7 | Medium | Pressure-level moist height uses actual rather than virtual temperature | Online/offline `:z` interpolation |
| F8 | Medium | Level-based PBL decay is pinned to 850 hPa | `:ModelLevel`; current `:PressureLevel` unaffected |
| F9 | Medium | Chunk rotation/restart can alter averaging intervals without metadata | Nonaligned cadence/restart; current cadence is aligned |
| F10 | Medium physical/numerical limitation | Thickness-weighted vertical momentum advection does not conserve kinetic energy | Stretched grids; current equal pressure thickness unaffected |
| F11 | Low for present primitive-equation work | Shallow-water output has incorrect physical units | Shallow-water analysis |

## F1 — Post-physics grid vorticity and divergence are stale

**Source:** [Spectral_Dynamics.jl](../../src/Dynamics/Spectral_Dynamics.jl),
`_synchronize_physics_next!`, lines 188–201; `Atmosphere_Update!`, lines
1092–1114; dynamics consumers at lines 695 and 801. See also
[Dyn_Data.jl](../../src/Core/Dyn_Data.jl), `Time_Advance!`, line 460.

The dynamics-only provisional state reconstructs `grid_vor` and `grid_div` from
`spe_vor_n` and `spe_div_n`. Physics then changes the winds, particularly through
Rayleigh friction. Synchronization recomputes the spectral vorticity/divergence
from these winds and reconstructs the winds, but never reconstructs their grid
vorticity/divergence. Time advancement copies the prognostic arrays and leaves
these diagnostic arrays unchanged.

Consequently, the next call uses post-physics winds together with pre-physics
divergence in mass continuity and pre-physics vorticity in the vector-invariant
momentum term. This is an equation-state inconsistency, not merely an output
problem. The `:vor` and `:div` output fields inherit it.

**Reproduction:** One normal primitive-equation step with Held–Suarez forcing
on T5L8 produces maximum errors, relative to inverse transforms of the accepted
current spectral state, of:

```text
vorticity: 9.465835046447e-9 s^-1
divergence: 1.431390215039e-10 s^-1
vorticity error / maximum current vorticity: 7.714868354892e-4
```

The same comparison with all physics disabled gives exactly zero vorticity
error. These numbers establish the defect; they are not an estimate of its
accumulated T42 climate impact. A small instantaneous discrepancy can still
modify repeated forcing and wave behavior over a long integration.

**Repair:** Reconstruct `grid_vor` and `grid_div` from the accepted post-physics
spectral fields after synchronization. Assert their consistency with the
spectral fields at every completed-step boundary in a short coupled regression.
Include a nonzero, horizontally varying wind subject to drag. The existing
coupling test verifies wind, temperature and pressure consistency but omits
these two fields.

## F2 — Two advertised tendency diagnostics contain only zeros

**Source:** [Spectral_Dynamics.jl](../../src/Dynamics/Spectral_Dynamics.jl),
`Reset_Dynamics_Tendencies!`, lines 112–125, and spectral tendency construction
around lines 806–835; [Variable_Mappings.jl](../../src/Core/Variable_Mappings.jl),
lines 128–129 and 150; [Output_Mappings.jl](../../src/Output/Output_Mappings.jl).

`grid_δvor` and `grid_δdiv` are allocated as zeros. Primitive-equation dynamics
resets them to zero, then calculates the actual tendencies in `spe_δvor` and
`spe_δdiv`. No source assignment transforms those tendencies to the grid
arrays. Output nevertheless maps `:dvor` and `:ddiv` to the zero arrays.
Barotropic `:dvor` has the same unpopulated diagnostic.

**Reproduction:** In the F1 step, maximum spectral vorticity and divergence
tendencies are respectively `3.936369380795e-10` and
`2.195818058529e-11 s^-2`; both grid tendency maxima are exactly zero.

**Impact:** These files silently suggest no vorticity/divergence tendency and
cannot be used for wave or process-budget diagnosis.

**Repair:** Define which tendency is intended—explicit dynamics, corrected
semi-implicit RHS, or actual completed-state change—and populate the grid
diagnostic accordingly. Document this choice. Existing generic `:du`, `:dv`,
`:dt`, `:dq`, and `:dps` working tendencies also do not collectively form a
complete physical tendency budget; their documentation already cautions against
that interpretation.

Two related diagnostic inconsistencies deserve regression checks:

- `grid_w_full` is calculated from the beginning-of-step state by
  `Four_In_One!` and is not rediagnosed after physics/time advancement.
  Saving it with accepted next-state winds/T/q introduces a one-base-step
  timing offset. This does not itself alter dynamics, which recalculates it.
- `grid_lnps` is reconstructed before the mass correction. That correction
  changes surface pressure and the mean spectral log-pressure coefficient,
  but does not refresh `grid_lnps`. `:lnps` output can therefore disagree with
  the accepted state. Refresh it or derive it from current pressure at output.

## F3 — A cold rerun can retain NetCDF chunks from the preceding run

**Source:** [Driver.jl](../../src/Driver.jl), cold-start cleanup at lines
293–303; [Output_Manager.jl](../../src/Output/Output_Manager.jl),
`_Init_Single_File`, line 109, and chunk construction/rotation.

Cold start removes old restart files but does not remove or reject an existing
NetCDF chunk series. The output manager replaces only the exact filename it
opens. A shorter rerun overwrites the early chunks it reaches while retaining
the old later chunks. Changing chunk frequency or pressure-output settings can
leave additional obsolete files as well.

**Reproduction:** An existing `output_t999999.nc` survives a fresh manager and
a complete shorter output sequence using the same base filename. The driver
has no NetCDF-series cleanup that would prevent this behavior.

**Impact:** A glob such as `output_t*.nc` can combine two physical simulations.
Contiguous timestamps alone need not detect the mixture. A “last 1000 days”
analysis may then read old results without an error. This is a potential
contamination path, not evidence that the already analyzed `ctrl`/`ctrl_BM`
directories actually contain mixed runs.

**Repair:** Prefer a new directory/run identifier for every cold experiment.
At cold start, reject conflicting existing chunks or archive the matching
series through an explicit workflow. Include a run ID and configuration/source
fingerprint in each chunk. Avoid broad automatic deletion of unrelated files.

## F4 — Restart loading is permissive about required state and configuration

**Source:** [Restart_Manager.jl](../../src/Output/Restart_Manager.jl), lines
41–49 and 69–81; [Driver.jl](../../src/Driver.jl), warm-start branch around
lines 260–282.

Version-3 files contain dimensions, time and array fields, but no timestep,
model type, vertical coefficients or planetary-constant fingerprint. Loading
checks dimensions only and skips every missing state field with
`haskey(file, path) || continue`. The driver immediately resumes mature
leapfrog with the newly configured timestep.

**Reproduction:** A version-3 file with valid dimensions/time and only
`state/grid_u_c` is accepted successfully. Required missing prognostics/history
retain their previous or newly allocated values. That loader success does not
establish a valid restart.

**Physical consequence:** If a checkpoint produced at timestep
`Δt_old` is resumed at `Δt_new`, its previous state is at
`t-Δt_old`, whereas leapfrog assumes it is at `t-Δt_new`. The update
`X_next = X_previous + 2Δt_new*RHS_current` has incompatible time history.
Changing the vertical coordinate with the same number of levels likewise loads
fields as if they occupied the new pressure surfaces, without remapping.
Changing model type can also pass the dimensions check.

**Repair:** Store and validate the integration timestep, equation set,
vertical coefficients and necessary physical constants. Define a required
prognostic/history schema while permitting missing optional diagnostics for
backward compatibility. If timestep changes are supported, explicitly rebuild
history/startup rather than silently continuing mature leapfrog. Intentional
physics-parameter perturbations need their own documented restart policy.

Existing restart tests pass array round-trip and legacy-format loading checks;
they do not establish trajectory equivalence across a resumed integration.
This finding concerns missing/incompatible state that is accepted.

## F5 — Several initialization paths do not create one consistent state

**Source:** [Initial_Conditions.jl](../../src/Initialization/Initial_Conditions.jl),
`Load_From_Arrays!`, lines 76–125, and `Init_Shallow_Water_Test!`, lines
244–246; [Variable_Mappings.jl](../../src/Core/Variable_Mappings.jl).

The array loader assigns `:ps` to `grid_ps_c`, but transforms `grid_lnps`
without first calculating `log(grid_ps_c)`. It assigns `:vor`/`:div` grid
fields without updating their spectral prognostics. When only one wind
component is supplied, the loop skips it and the paired-wind branch never
executes, silently ignoring the supplied component. Unsupported keys are
also silently skipped.

**Reproduction:** Supplying uniform `:ps = 100000 Pa` to a fresh
primitive-equation state leaves the corresponding spectral pressure equal to
`exp(0) = 1 Pa`. Supplying nonzero barotropic `:vor` leaves `spe_vor_c` zero.
The array loader also does not initialize previous-level conservation/history
arrays. That history limitation is acknowledged in `docs/src/config.md`,
but it compounds the incorrect pressure alias.

There are two additional 2-D problems:

- The configuration uses `:Shallow_Water`, but the variable map dispatch is
  `Val(:ShallowWater)`. Array initialization of the configured shallow-water
  model throws a `MethodError` before loading. This part is loud rather than
  silent.
- Built-in shallow-water initialization creates uniform `grid_lnps` and
  `spe_lnps_c = h_0`, but leaves the actual `grid_ps_c` height field zero.
  For `h_0=30000`, the first physics evaluation therefore relaxes from zero
  grid height while dynamics uses spectral height 30000. The probe confirms
  those exact values. Built-in barotropic initialization similarly adds a
  vorticity perturbation without reconstructing the corresponding initial
  winds, so its first nonlinear tendency uses the unperturbed winds.

**Repair:** Centralize initialization finalization: select authoritative
prognostics, perform pressure/log-pressure and height alias conversion,
project/reconstruct all corresponding grid and spectral fields, and initialize
required history. Normalize the model-type symbol and reject unknown or
incomplete array specifications. These findings do not apply to the present
`Moist_Spinup` initialization, whose grid/spectral/current/previous fields are
explicitly initialized and tested.

## F6 — The energy-correction denominator omits finite model-top pressure

**Source:** [Spectral_Dynamics.jl](../../src/Dynamics/Spectral_Dynamics.jl),
`Compute_Corrections!`, lines 324–342; mass definition in
[Vert_Coordinate.jl](../../src/Core/Vert_Coordinate.jl).

The energy integral uses actual layer pressure masses, but a uniform
temperature correction is calculated as

```text
ΔT = g (E_target-E_next) / (cp * mean(ps)).
```

For a fixed finite top, the represented atmospheric mass is
`mean(ps-p_top)/g`. The correct denominator for this *same defined ledger*
is `cp * mean(sum(Δp))`, not `cp * mean(ps)`.

**Reproduction:** A 100000 Pa surface, 10000 Pa fixed top, and uniform 1 K
temperature perturbation retain **0.1 K**, or **10% of the energy discrepancy**,
after the purported correction. The result agrees with the exact missing-mass
factor `p_top/ps`. The supported Simmons–Burridge table has a 219.422 Pa top;
the default hybrid grid has a finite top near 1831.56 Pa. The residual is
smaller on these grids but still nonzero.

**Repair:** Use the global represented column mass derived from the same
`Δak + Δbk*ps` used by the integral. Test the correction independently on a
finite-top grid. This repairs enforcement of the existing energy ledger; it
does not prove that the ledger includes every physical boundary-energy term.
Current zero-top even sigma does not activate this denominator error.

## F7 — Moist pressure-level height interpolation uses a dry hypsometric term

**Source:** [Vertical_Interpolation.jl](../../src/Output/Vertical_Interpolation.jl),
lines 269–277; [Output_Manager.jl](../../src/Output/Output_Manager.jl),
lines 449–466; [Interpolator.jl](../../post_processing/Interpolator.jl),
temperature loading and interpolation around lines 193–218.

Native heights are integrated hydrostatically using virtual temperature:

```text
Tv = T * (1 + (Rv/Rd-1)*q).
```

Pressure-level interpolation instead adds a hypsometric height difference
using the actual-temperature array. Neither caller supplies virtual
temperature or humidity to that calculation. Moist native and pressure-level
height therefore use different equations of state, even for an isothermal
column. The recently added default offline `:z` output exercises this path.

**Reproduction:** For constant `T=300 K`, `q=0.02`, native nodes at 500/900
hPa and a target at 700 hPa, the interpolated height is **26.84346 m too low**
relative to the analytic moist hydrostatic height. Dry humidity eliminates
the mismatch.

**Repair:** Supply an equation-of-state temperature to height interpolation
and carry the required moisture information through online/offline output.
For averaged records, explicitly define how the virtual-temperature average
is obtained: `mean(Tv)` includes a `mean(T*q)` term and is generally not equal
to `Tv(mean(T),mean(q))`. Ordinary log-pressure interpolation of other variables
is a separate operation. This finding changes diagnostics, not native dynamics.

## F8 — Model-level PBL top does not control the decay-profile pressure

**Source:** [PBL.jl](../../src/Physics/PBL.jl), lines 452–463.

The `:ModelLevel` option determines the interfaces with constant mixing, but
the exponential decay above them uses a hard-coded 85000 Pa rather than the
pressure of the selected top interface. The effective decay can therefore be
discontinuous and almost absent just above a deep selected PBL.

**Reproduction:** On a 20-level even-sigma column with eight bottom layers
selected, the constant-mixing top is at 600 hPa. At the next interface, 550
hPa, the code gives

```text
K(550)/K_base = exp(-((850-550)/100)^2) = 0.0001234098.
```

A decay anchored at the selected 600 hPa top gives `exp(-0.25) = 0.7788008`.
The discrepancy is approximately a factor of 6300.

**Repair:** Obtain the tail's top pressure from the selected column interface,
and validate the selected level range. The current pressure-based 850 hPa
configuration uses the correct branch and is unaffected by this finding.

## F9 — Output intervals can change silently at chunk/restart boundaries

**Source:** [Driver.jl](../../src/Driver.jl), cadence validation at lines
104–113; [Output_Manager.jl](../../src/Output/Output_Manager.jl), lines
331–345, 388–407, and 543–576; restart format described in F4.

Fields are accumulated every timestep and averaged over `sample_counter`.
They are **time means**, with an ending-time coordinate. NetCDF variables lack
`cell_methods="time: mean"`, and no time bounds or record sample count are
written. Under CF conventions, means should be identified and their intervals
described; the files declare CF-1.11. See
[CF cell methods](https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/ch07s03.html).

Validation requires output and saving cadences to be divisible by `Δt`, but
does not require chunk boundaries to align with output intervals. Rotation
flushes a partial mean, and ordinary flushes remain aligned to the original
start time. Finalization also flushes a partial mean. Restart saves no output
accumulator and starts a new averaging schedule at the checkpoint time.

**Reproduction:** With 600 s timesteps, a 3600 s output interval and a chunk
boundary at 5400 s, record ending times are `3600, 5400, 7200, 10800 s`.
The middle two records cover 1800 s each while the others cover 3600 s. For
the test field `u(time)=time`, recorded means are respectively
`2100, 4800, 6600, 9300`, confirming the different sample windows.

**Impact:** A record-count average weights unequal intervals incorrectly.
Uniform-cadence FFT assumptions fail, and files do not contain enough metadata
to reconstruct all averaging windows reliably. Restarting within a window
also changes the result relative to uninterrupted output.

**Repair:** Preserve windows across chunk boundaries, or require aligned
cadences and mark any partial final window. Write time bounds and
`cell_methods`, and save accumulator state when exact output continuation is
required. The current 50-day/12-hour cadence is aligned, so the rotation defect
does not automatically contaminate those records. Their interpretation as
12-hour means remains essential. My earlier description of them as snapshots
was incorrect.

## F10 — Weighted vertical momentum interpolation violates the kinetic-energy identity

**Source:** [Vert_Coordinate.jl](../../src/Core/Vert_Coordinate.jl), lines
618–641 and 658.

This is a confirmed discrete conservation limitation on stretched grids. For
a closed vertical column, let `m_k = Δp_k`, with continuity
`m_dot = M_top-M_bottom`. A consistent centered momentum transport should
satisfy the vertical-transport kinetic-energy identity

```text
sum(m_k*u_k*u_dot_k + 0.5*u_k^2*m_dot_k) = 0.
```

At an interior face, the contribution is
`M*(u_lower-u_upper)*(u_face-(u_upper+u_lower)/2)`.
Arithmetic averaging makes it zero even when adjacent layer masses differ.
The thickness-weighted face value does not. The sign reverses with the flux,
so it can numerically create as well as remove kinetic energy.

**Reproduction:** Two layers with pressure thicknesses 10000/30000 Pa,
winds 10/20 m/s, a 1 Pa/s downward interior flux, and closed boundaries give
a pressure-weighted kinetic-energy rate of **-25 Pa m² s^-3** with
`second_centered_wts`, versus **zero** with `second_centered`. Divide by
gravity for the physical rate per unit area.

**Impact:** The global energy fixer can hide this exchange by applying a
uniform temperature correction. That keeps the chosen global ledger close
while altering the spatial distribution of momentum/heating. Spatial
interpolation accuracy and energy conservation are distinct requirements;
this finding does not mean that every weighted interpolation is unusable.

**Repair:** Decide and document the intended discrete momentum/energy
conservation property. Use the arithmetic face value or a compatible
energy-conserving formulation, then test its local energy identity independently
of the global fixer. Equal-thickness even sigma makes the two face averages
coincide, so this defect is not active in the present T42 vertical grid.

## F11 — Shallow-water variables are geopotential but labeled as height

**Source:** [Shallow_Water_Dynamics.jl](../../src/Dynamics/Shallow_Water_Dynamics.jl),
energy construction around line 160; [Output_Mappings.jl](../../src/Output/Output_Mappings.jl),
lines 184 and 190–201; [Variable_Mappings.jl](../../src/Core/Variable_Mappings.jl),
lines 158–161.

The shallow-water scalar is added directly to kinetic energy and differentiated
as a pressure-force potential. Its units are `m² s^-2`, but `:h` is exported
in `m`, and `:dh` in `m s^-1`. The existing dynamics documentation acknowledges
this unresolved metadata mismatch. Potential vorticity is divided by this
geopotential scalar while its metadata assumes division by height.

**Repair:** Either retain geopotential and give it appropriate names/units, or
divide by gravity when exporting height, height tendency and height-based PV.
PV also uses `grid_absvor`, computed before advancement, with next-state height;
calculate current absolute vorticity when diagnosing current PV.
This does not affect the primitive-equation WK results.

## Physics and numerical checks without a newly established defect

- Gaussian quadrature, Legendre normalization, Laplacian signs, wind
  reconstruction and spectral transform normalization are consistent with the
  tested spherical-harmonic identities. The configuration grid constraint
  provides the intended quadratic-product sampling for ordinary truncations.
- The semi-implicit Helmholtz inverse has the expected sign with negative
  Laplacian eigenvalues. The back-substitution and mature-step matrix update
  follow the stated leapfrog equations. Its dry 300 K linear reference is an
  approximation to moist dynamics, not evidence of an algebraic sign error.
- The saturation formulas distinguish mixing ratio from specific humidity;
  their derivatives include the phase blend. Virtual temperature is used in
  native moist hydrostatics, pressure gradients and pressure work.
- BM ascent, LCL handling, virtual buoyancy, RH conversion, deep precipitation
  closure and fractional shallow adjustment have the prior independent Isca
  comparisons and passing column tests. No new ordinary-path column-algebra
  defect was found. Edge-domain caveats remain as described in
  [the earlier BM audit](static_audit_20261005.md).
- Condensation is a documented linearized saturation adjustment, and can
  differ from solving the nonlinear saturation constraint exactly. It is not
  an exact saturation solver. Its finite increments and precipitation signs
  are consistent with the implemented formulation.
- Physics covers one base interval at startup and two bounded substeps during
  mature leapfrog. A `2Δt` physics interval is therefore not, by itself, double
  application: the provisional dynamics state also connects previous and next
  levels over `2Δt`.
- Rayleigh friction returns its removed kinetic energy to temperature. PBL
  mixes potential temperature and humidity; low-level momentum drag supplies
  the momentum representation. Absence of a separate PBL wind-diffusion solve
  is consistent with the idealized moist Held–Suarez design described by
  [NCAR](https://www.cesm.ucar.edu/models/simple/moist-held-suarez).
- Grid moisture transport is in advective form, with fixed winds/masses per
  substep, horizontal/vertical splitting and positivity limiting. A global
  water correction is used to reconcile its integral with the core; the raw
  transport should not be described as an exactly locally mass-coupled moist
  continuity solver. The constant-field and high-Courant positivity checks pass.
- The regularized LRF contracts, vertical matrix indexing, degree-valued
  latitude comparison and K/day-to-K/s conversion are consistent with the
  current tests. Its loader checks dimensions and latitudes but has no model
  vertical-coordinate argument. Artifact/grid pressure correspondence remains
  an external requirement, especially if the number of levels stays unchanged
  while their pressure locations change.
- The filter is ordinary Robert–Asselin, despite residual RAW/Williams wording
  in source comments. `docs/src/dynamics.md` already states this accurately.
  Switching to RAW would be a numerical-method change, not a one-line
  correction. The distinction and physical-mode damping are explained by
  [Williams (2009)](https://journals.ametsoc.org/view/journals/mwre/137/8/2009mwr2724.1.xml).

## Physical limitations relevant to the ±5 signal

The code imposes zero relative mass flux at the upper boundary and has no
dedicated upper sponge. Held–Suarez momentum drag operates near the surface;
horizontal hyperdiffusion is not an upper absorbing layer. Reflection and
trapping remain hypotheses to test, rather than proof of an implementation
error. Twenty uniformly spaced sigma layers also provide limited vertical
resolution aloft.

Sudden convective/condensation adjustments and sequential dynamics–physics
coupling can excite fast gravity waves even with a correct column scheme.
This phenomenon is described in
[Thatcher and Jablonowski (2016)](https://gmd.copernicus.org/articles/9/1263/2016/).
F1 introduces an additional, independently demonstrated state inconsistency;
both causes can coexist.

Global energy/water corrections enforce selected scalar totals and can conceal
local numerical errors. The earlier five-day ledger groups dynamics and its
fixers, and explicitly measures nonzero PBL/pressure-adjustment energy effects.
Its small total residual does not test the F1 state invariant or establish a
complete moist-energy theorem. The default prescribed SST and Newtonian
forcing are also idealized boundary/thermal choices, not full radiative or
ocean energy closure.

The stored 12-hour means have a Nyquist frequency of 1 cycle/day. Averaging
attenuates high frequencies but does not independently eliminate aliasing.
Short periods close to this limit require higher-cadence output for confident
wave identification. None of the static findings uniquely selects zonal
wavenumber 5.

## Recommended repair and validation order

1. Repair F1 and add a completed-step grid/spectral vorticity/divergence
   invariant. Diagnose and record pre-fixer mass, energy and moisture changes
   in a short coupled run so global corrections cannot conceal the result.
2. Repair the zero and lagged diagnostics before using tendency/phase budgets.
   Separate accepted-state diagnostics from interval/process-rate diagnostics.
3. Address run-series integrity and restart compatibility. Use fresh output
   directories for every subsequent paired experiment.
4. Repair and independently test finite-top energy correction, initialization
   aliases, moist height interpolation and model-level PBL selection. Establish
   the intended momentum-energy discretization before adopting stretched grids.
5. Rerun matched BM-off/BM-on experiments with identical current sources,
   initial states and sampling windows. Save hourly or finer wind, temperature,
   humidity, pressure velocity and named convective tendencies for a short wave
   diagnosis. Compare raw power as well as normalized WK power.
6. Then vary upper damping and vertical resolution, and separately halve the
   timestep/modify physics application cadence. These controlled experiments
   distinguish reflected waves, coupling excitation and physical moist modes.

A comprehensive static audit cannot certify a long-run climate or attribute a
specific spectrum peak. The concrete outcome here is an independently
reproduced active coupling defect, several conditional implementation/output
defects, and a bounded list of physical/numerical limitations to test.
