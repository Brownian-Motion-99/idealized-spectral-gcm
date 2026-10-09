# LRF implementation and experiment readiness, 2026-10-08

This records the October 8 implementation and experiment snapshot. The current
[physics documentation](../../docs/src/physics.md#moisture-linear-response-function)
now covers all three artifact schemes. The main HSt42 launcher also uses an
LRF-specific output suffix and an explicit experiment SST callback.

## Assessment

The current working-tree implementation has a complete path from coefficient
loading to temperature increments and NetCDF diagnostics. The
latitude-specific logarithmic scheme is the best supported candidate for
controlled T42L20 experiments. This is a conditional assessment: the supplied
artifact is calibrated against the older `ctrl` climate, the experiment script
needs separate output paths and an intentional initialization choice, and short
numerical checks cannot establish equilibrium climate stability or accuracy
after a large change in climate.

In particular, the current artifact adds approximately **+5.2 W/m²** global
column heating in the first-day `ctrl_BM` pilot. It is therefore not established
as a feedback centered on the BM control. That experiment objective requires
recalibration/validation; an experiment deliberately retaining the original
common reference can instead treat this mean adjustment as part of its forcing.

The shared tropical matrix with the supplied 26°/45° taper is less well
supported physically. Its substantial errors at 15–26° are demonstrated in an
independent cached validation. The legacy linear-humidity reader remains
available, with weaker input validation and no new physical certification.

This review covers the **current working tree**, including its existing
uncommitted LRF and post-physics synchronization changes. It does not describe
only the committed version. Production source, existing experiments, external
artifacts, and control outputs were left unchanged. The additions are this
report, review evidence scripts, and small recorded results.

## What the feature computes

LRF is an additive, prescribed longwave temperature response to humidity. It
changes temperature directly; it does not directly change humidity, winds, or
surface pressure. Moist physics and dynamics can respond to the temperature
change in later substeps and steps. It is used alongside Held–Suarez Newtonian
relaxation in the provided experiments.

There are three runtime states in [LRF.jl](../../src/Physics/LRF.jl):

| Artifact scheme / Julia type | Heating before conversion to K/s | Required arrays |
| --- | --- | --- |
| `linear_q_v1` / `LRF_State` | `L(latitude) * (q - q_reference)` | `LRF_LW_q(nd,nd,nlat)`, `ref_q(nlon,nlat,nd)` |
| `regularized_log_tapered_v1` / `Regularized_LRF_State` | `alpha * taper(latitude) * B_star * (log(q+q0) - chi_reference)` | `B_star(nd,nd)`, `chi_reference(nd,nlat)`, `taper(nlat)` |
| `regularized_log_latitude_v1` / `Latitude_LRF_State` | `alpha * B_by_lat(:,:,latitude) * (log(q+q0) - chi_reference)` | `B_by_lat(nd,nd,nlat)`, `chi_reference(nd,nlat)` |

For both logarithmic schemes the file also supplies `q0_kg_kg`, `alpha`, and,
on the driver path, `latitude` in degrees. A file without `scheme` is interpreted
as the legacy scheme. These are JLD2 runtime artifacts; the offline HDF5 files
are not drop-in runtime inputs.

Specific humidity is in kg/kg. The legacy matrix is interpreted as K/day per
unit specific humidity. The logarithmic matrices are K/day per unit natural-log
humidity perturbation. All schemes divide their result by `config.day_to_sec`
to return K/s; the current experiments use 86,400 s/day. Matrix rows identify
the temperature response level and columns identify the humidity perturbation
level, so the response can couple all levels within a column. The supplied
artifacts are ordered from the atmosphere's top to the surface.

The kernel is evaluated independently at every longitude and latitude: there is
no horizontal convolution and no zonal averaging of the online humidity. Each
latitude-specific matrix is shared by all longitudes at that latitude. The
legacy reference may depend on longitude; logarithmic references depend on
level and latitude only. There is no online latitude or vertical interpolation.

The logarithmic variants are linear in the transformed humidity anomaly, not
in specific humidity itself. Locally, the humidity derivative of each matrix
column scales as `1/(q+q0)`, together with `alpha` and any taper. The positive
`q0` makes exact zero humidity finite; it is not a humidity floor applied to the
prognostic state or a cap on heating. Negative or nonfinite humidity is rejected.

The supplied files use `q0=1e-8 kg/kg` and `alpha=1`. The tapered file has full
strength through |latitude|=26°, a cosine transition, and zero feedback at
|latitude|>=45°. These choices are stored in the artifact; the runtime consumes
the stored taper rather than generating it. The latitude-specific file has a
20×20 matrix for each of the 64 Gaussian latitudes and no prescribed taper.

## Calibration and fixed reference

The external builder
[/home/garywu/undergrad_proposal/LRF/build_latitude_lrf.py](/home/garywu/undergrad_proposal/LRF/build_latitude_lrf.py)
calculates a humidity Jacobian of clear-sky RRTMG longwave heating at each
latitude's control mean column. Atmospheric temperature, SST, nominal pressure,
and other absorbers remain fixed while humidity is perturbed. The selected
centered finite-difference step is ±10%; ±5% and ±20% are convergence checks.
The transformed matrix scales the specific-humidity Jacobian columns by
`q_calibration+q0`.

For these files, `chi_reference` is the time/zonal mean of `log(q+q0)` in
`/data92/garywu/undergrad_proposal/ctrl`, days 2650.5–3650, with 2,000 archived
times. It is not `log(mean(q)+q0)`. Its equivalent humidity reference is
`exp(chi_reference)-q0`. This centering makes the time/zonal mean added heating
zero in the calibration sample at each latitude and level, to roundoff.

Matrices, reference, strength, and taper stay fixed throughout integration.
The feature does not call RRTMG online, evolve a radiation state, predict cloud
or shortwave feedback, or include a temperature/SST radiation response in the
LRF itself. It adds humidity-induced anomalous heating rather than the stored
reference RRTMG heating. The latter is not added to Held–Suarez forcing.

The artifact remains tied to nominal reference pressure rather than changing
column pressures. These are limitations of the prescribed experiment design,
not matrix indexing defects.

## Loading, coupling, and diagnostics

[Driver.jl](../../src/Driver.jl) copies `physics_params` and loads the file once
before time integration when `"do_LRF"=>true`. It requires moisture transport
and a string `"LRF_file"`. The loaded state is placed in the runtime dictionary
under `"LRF_state"`; no coefficient file is read in the timestep loop.

The loader checks required keys and array dimensions. Both logarithmic
constructors reject nonfinite matrices/references, invalid `q0`, negative or
nonfinite `alpha`, and, for the tapered variant, weights outside [0,1]. On the
driver path their stored latitude values and ordering must match
`rad2deg.(mesh.θc)` to an absolute tolerance of 1e-10 degrees. Direct loader
callers may omit this coordinate check by leaving `latitude=nothing`.

[Spectral_Physics_Interface.jl](../../src/Physics/Spectral_Physics_Interface.jl)
applies physics to the provisional next state from dynamics. Each base substep
uses this order:

1. Betts–Miller convection and large-scale condensation;
2. surface sensible heating and evaporation;
3. PBL mixing;
4. Rayleigh friction and Newtonian relaxation;
5. LRF heating;
6. physical-state validation and dry-air/pressure adjustment.

Thus LRF sees humidity after the moisture processes and PBL mixing, before the
end-of-substep dry-air adjustment. Its explicit increment is
`T += physics_dt * lrf_temperature_rate`. Startup has one 600 s substep in the
current configuration. Mature leapfrog covers the 1,200 s previous-to-next
interval with two 600 s physics substeps. The LRF diagnostic is the mean of the
substep rates. This is not a second insertion into the leapfrog RHS.

After physics, [Spectral_Dynamics.jl](../../src/Dynamics/Spectral_Dynamics.jl)
projects the adjusted state into the spectral representation and applies global
corrections using the **post-physics** energy target. The energy fixer therefore
does not impose the pre-LRF energy integral and cancel the physical heating.
The current working tree also refreshes post-physics grid vorticity/divergence;
the stale-field defect described in the earlier October 6 audit has been fixed.

`LRF!` overwrites its dedicated output buffer and leaves its humidity input
unchanged. Logarithmic states keep one reusable anomaly column per Julia
thread, and kernels thread over latitude. Arithmetic cost is
`O(nlon*nlat*nd^2)` per evaluation. State and diagnostic arrays persist between
steps; there is no per-column full-grid allocation or runtime radiation call.
The static threading and shared scratch storage assume the normal sequential
driver invocation rather than concurrent calls on the same state.

The live diagnostic is `dyn_data.grid_lrf_tendency`. Requesting `:lrf_dt`
writes `lrf_dta_dt`, a 3-D temperature tendency in K/s. Diagnostics are reset
each dynamics call, including when LRF is disabled. Output records average
these rates over the configured output interval. They are process tendencies
before spectral projection/correction and should not be equated point by point
with the entire accepted-state temperature difference.

Restart files save model arrays and time. LRF coefficients are loaded again
from the externally selected file when restarting; they are not embedded in the
checkpoint. NetCDF global attributes currently do not identify the LRF scheme,
artifact hash, reference period, `alpha`, or `q0`.

## Findings that affect experiments

### 1. The main T42 launcher reuses the control output path

`exp/HSt42/HS.jl` enables LRF whenever `JGCM_LRF_FILE` is nonempty, but still
sets `experiment_name="ctrl_BM"`, uses the same control output directory, and
sets `is_restart=false`. The cold-start driver removes JLD2 files in that
directory's `restart` subdirectory. Existing NetCDF chunks can also be
overwritten or survive alongside a new shorter run. Simply supplying the
environment variable to this unchanged launcher is unsuitable for retaining
a clean control/treatment pair.

Use a distinct experiment/output directory and a deliberately selected warm
restart. This is a required experiment configuration change, not an LRF kernel
repair. The review integrations use fresh directories under `/tmp` and never
write to the control directories.

### 2. Vertical pressure and artifact provenance are unchecked

The loader has no model vertical-coordinate argument and ignores the stored
`pfull`, `phalf`, source hash, and cross-language validation arrays. A file
with reversed pressure metadata is accepted if its array dimensions and
latitudes match. A different 20-level grid can therefore be silently paired
with these 20-level matrices.

Both supplied artifacts match the current even-sigma nominal pressure grid:
maximum full-level discrepancy is 1.14e-13 hPa; interface discrepancy is zero.
The missing guard does not invalidate this verified configuration. Changing
vertical grid, ordering, units, or artifacts requires a new explicit check.

### 3. The tapered shared matrix has a demonstrated regional accuracy limit

The cached independent 960-column validation at |latitude|=15.35–23.72° gives
the exact supplied tapered formulation an observed-anomaly median profile
error of **61.5%**, with a **160.4%** 90th percentile. Lower-tropospheric
moistening/drying median errors are 75.7%/86.0%. This region is inside its
full-strength 26° core, so the stored taper does not reduce those errors.
Numerical tests passing does not support treating this artifact as an accurate
response throughout its active latitude range.

The latitude-specific artifact gives substantially better cached held-out
performance over 1,536 sampled control columns:

| Perturbation | Median relative profile error | 90th percentile | Median profile cosine |
| --- | ---: | ---: | ---: |
| Lower-level moistening | 16.2% | 37.9% | 0.992 |
| Lower-level drying | 14.9% | 38.6% | 0.993 |
| Free-tropospheric moistening | 20.3% | 39.1% | 0.985 |
| Free-tropospheric drying | 20.9% | 37.4% | 0.984 |
| Observed humidity anomaly | 15.0% | 26.8% | 0.993 |

These statistics were recomputed from the cached HDF5 metrics during this
review; new RRTMG calculations were not run. The independent window is days
2000.5–2500. The maximum matrix step sensitivity is 1.081%, below the 1.5%
configured convergence tolerance. Held-out mean-heating RMS is 0.0311 K/day,
global signed column heating is -0.0837 W/m², and disturbance RMS is
0.5323 K/day. These comparisons hold temperature/SST fixed and do not validate
equilibrated climate response or transfer to a warmer climate.

### 4. The supplied reference belongs to `ctrl`, not a demonstrated `ctrl_BM` calibration

The current main experiment turns on Betts–Miller with RH=0.7 and tau=7,200 s.
The artifact records the older `ctrl` baseline. It does not establish zero
time-mean added heating around the new `ctrl_BM` climatology or a changed SST
climate. A mean adjustment around those states may be an intended consequence
of using a common fixed reference; it should be measured and stated in the
experiment design.

The paired one-day pilot makes this concern concrete: LRF's global column
heating in successive six-hour records is +5.248, +5.214, +5.183, and
+5.126 W/m², using nominal layer mass and Gaussian area weights. Its one-day
time/zonal-mean heating profile has 0.408 K/day mass/area-weighted RMS. These
short-window values establish an appreciable mean adjustment during the pilot;
they do not estimate a stationary climatological bias.

If the intended treatment is centered on `ctrl_BM`, first estimate its mean
`log(q+q0)` and validate the humidity response on that control. If the intention
is to retain the original common reference across experiments, preserve it and
diagnose the resulting mean forcing. Automatically recentering every treatment
would change the experiment's definition.

### 5. Zero-safe logarithms do not ensure weak cold-start forcing

There is no online activation ramp, heating limiter, or automatic recalibration.
Evaluating the supplied latitude-specific artifact on the current analytical
cold-start humidity gives rates from -1.093 to +3.580 K/day before other physics.
The tapered artifact gives -1.489 to +1.765 K/day. These current values should
not be replaced by the older README's different cold-start extreme.

An all-zero humidity stress input remains finite but yields latitude-specific
rates from -17.45 to +27.04 K/day. This demonstrates what regularization does
and does not guarantee; it is not a representative equilibrium climate.
Using an equilibrated matching control restart is the appropriate default
for an experiment applying these fixed reference anomalies.

### 6. Legacy validation and experiment controls are weaker

`LRF_State` accepts nonfinite coefficients/reference values; its kernel can
produce NaN before the general physical-state validator stops integration.
The legacy loader also skips latitude comparison even when model latitudes
are supplied. Both behaviors were reproduced in isolated temporary fixtures.

For logarithmic schemes, `alpha` and `q0` are read from the artifact. There is
no implemented `physics_params["LRF_alpha"]` override or activation schedule.
A strength sweep currently requires separately recorded artifact variants or
deliberate state construction outside the standard file-loading driver path.
Keep the selected runtime file and its hash with each experiment because the
restart/output files do not establish which external coefficients were used.

At review time the public physics documentation described only the legacy
formulation and file keys. The current
[physics documentation](../../docs/src/physics.md#moisture-linear-response-function)
includes the logarithmic formulations, artifact layouts, and coupling.

## Verification and reproduction

Julia 1.8.5 was used. The existing complete regression suite passes **16,651**
assertions: 90 core numerics, 170 runtime infrastructure, 15,045 physics, and
1,346 dynamics–physics integration assertions. The dedicated LRF file passes
all **51** assertions with four threads.

The new evidence harness passes **36** assertions with four threads for the
actual artifacts, reference heating, zero humidity, alpha=0, unchanged humidity
input, Python/Julia agreement, nominal pressure correspondence, both physics
intervals, and preservation of the post-physics energy target. Some assertions
explicitly demonstrate accepted invalid metadata and legacy validation gaps;
they are evidence of limitations rather than requirements for future behavior.

Maximum Python/Julia heating errors are 2.75e-20 K/s for the latitude-specific
file and 2.77e-20 K/s for the tapered file. Both match their embedded validation
column to roundoff.

The full driver completed **1,008 T42L20 steps** across three warm-start runs,
with 600 s base timestep and **30 saved-output assertions** passing:

| Run | Duration and control checkpoint | Saved-grid LRF RMS, K/day | Saved temperature range, K |
| --- | --- | ---: | ---: |
| Latitude LRF, BM off | 5 days from `ctrl`, day 3650 | 0.368–0.394 | 168.85–300.42 |
| LRF off, BM on | 1 day from `ctrl_BM`, day 3650 | Exactly zero | 183.86–299.83 |
| Latitude LRF, BM on | 1 day from the same `ctrl_BM` checkpoint | 0.497–0.512 | 184.39–299.76 |

All saved winds, temperature, humidity, surface pressure, and LRF diagnostics
are finite. Saved humidity is nonnegative and below one; surface pressure is
positive. The first run saves daily means; the paired BM runs save six-hour
means. These are checks of saved interval means, not of instantaneous extreme
heating or long-term climate equilibration. The normal model state validators
also run during integration.

The paired runs differ in their final six-hour temperature means by 0.229 K
grid RMS, with maximum absolute difference 2.584 K. These differences include
the dynamical/moist-physics response and are not a direct LRF tendency budget
or an equilibrium sensitivity estimate.

Recorded evidence:

- [Actual-artifact contracts and source hashes](lrf_review_20261008_contracts.toml).
- [Three online runs and output bounds](lrf_review_20261008_online.toml).
- [Paired comparison and approximate mean column heating](lrf_review_20261008_comparison.json).
- Regression-suite counts are recorded above; the suite can be rerun below.
- [Julia evidence harness](lrf_review_20261008.jl) and
  [Python comparison utility](lrf_review_20261008_compare.py).

```sh
julia --project=. --compiled-modules=no --threads=1 test/runtests.jl
julia --project=. --compiled-modules=no --threads=4 \
  -e 'using Test, JGCM; include("test/test_LRF.jl")'
julia --project=. --compiled-modules=no --threads=4 \
  test/validation/lrf_review_20261008.jl contracts /tmp/NEW_LRF_CONTRACTS
julia --project=. --compiled-modules=no --threads=1 \
  test/validation/lrf_review_20261008.jl online /tmp/NEW_LRF_ONLINE
```

The harness requires a new output directory and reads the existing external
artifacts and day-3650 `ctrl`/`ctrl_BM` checkpoints. Machine-readable results
contain relevant production source hashes and artifact hashes. Large generated
NetCDF data stay under `/tmp`; the production 3,650-day integrations were not
launched.

The five-day model run was completed under `/tmp/lrf_review_20261008_online/`.
An initial review-reader failure converting CF-decoded dates to Float64 was
corrected to read numeric time values; the completed NetCDF was then reused
with `JGCM_REVIEW_BM_OFF_NC` rather than repeating the model run. Both BM pilots
and final online results are under `/tmp/lrf_review_20261008_online_v2/`.
The separate comparison was made with:

```sh
/home/garywu/.conda/envs/climlab/bin/python \
  test/validation/lrf_review_20261008_compare.py \
  /tmp/lrf_review_20261008_online_v2/control_bm_on_1day/output_t315360000.nc \
  /tmp/lrf_review_20261008_online_v2/latitude_bm_on_1day/output_t315360000.nc \
  test/validation/lrf_review_20261008_comparison.json
```
