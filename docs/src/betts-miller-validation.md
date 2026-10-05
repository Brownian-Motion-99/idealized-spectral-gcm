# Betts--Miller validation

The scientific target is Frierson's Simple Betts--Miller formulation with the
O'Gorman--Schneider modifications in Isca's `qe_moist_convection.F90`.
Validation separates the fixed-mass column adjustment from its full-model
coupling and integration behavior.

## Static audit and independent columns

The 2026-10-05 static audit of implementation `e0ba262` found no further error
in the ordinary deep or shallow adjustment algebra. It checked parcel ascent,
LCL bounds, virtual buoyancy, exact reference RH, both deep closures,
fractional shallow adjustment, units, scratch storage, and process coupling.
The audit records CIN's signed convention, compatibility handling of negative
standalone humidity, and the limits of the column energy identity.

The local Isca reference revision is
`5c6d5baf74498c5dce3578e235c5254e6856baf6`. An isolated double-precision
Fortran harness matches constants, the direct saturation kernel, virtual
temperature, RH 0.8, tau 7200 s, and minimum parcel temperature 173.16 K.
Of 114 validation inputs, 87 are compared with the original Isca code; dry
and model-top edge cases use independent physical/indexing invariants.
The comparison covers both deep branches and shallow transport.

The largest temperature-rate difference from the original Isca LCL lookup is
`1.77e-7 K/s`. A separate control replaces only that lookup with log-pressure
bisection, reducing the difference to `1.96e-17 K/s`. Original interpolation,
surface-coincident LCL rounding, and source-reference boundary residues are
explicitly reported, rather than treated as exact equivalence.

The BM-focused suite passes 15,024 assertions on one and four threads. Tests
cover water/rain and moist-enthalpy budgets, finite steps through `dt=tau`,
mixed scalar/grid agreement, work-buffer reuse, and startup/leapfrog coupling.
Full-suite validation passes 15,626 assertions when the pre-existing default
SST edit is excluded in a disposable copy. That edit causes three existing
prescribed-temperature failures in the original working tree and was preserved.

The detailed audit, fixed Fortran fixtures, optional comparison utility, and
source hashes are in `test/validation/` and `test/fixtures/` in the repository.
These normal tests need neither Isca nor a Fortran compiler.

### Before and after the corrections

The same fixed inputs were evaluated with the original Julia implementation
at `6197034`, in an isolated archive, and the corrected implementation at
`e0ba262`. Constants and the virtual-temperature setting are identical; the
old implementation ignores that setting in its buoyancy calculation.

| Diagnostic | Original | Corrected |
|:---|---:|---:|
| Near-neutral moist column CAPE, J/kg | 0 | 454.23657 |
| Shallow column CAPE, J/kg | 586.25894 | 988.82923 |
| Shallow maximum absolute humidity rate, s⁻¹ | 0 | 4.78589e-7 |
| Shallow rain, mm/day | 0 | 0 |
| Deep column rain, mm/day | 44.64631 | 48.16968 |
| Deep surface temperature rate, K/s | approximately 0 | -1.47822e-4 |
| Dry-neutral maximum parcel temperature error, K | 0.243738 | 0 |
| Dry-neutral CAPE, J/kg | 0.547973 | 0 |

The deep-rain change combines buoyancy and RH corrections and should not be
attributed only to the closure change. The shallow case now transports heat
and moisture with zero rain. The original RH multiplier gives actual RH
`0.80569062` at 300 K and 1000 hPa for requested RH 0.8; the corrected
vapor-pressure target gives 0.8 to roundoff.

Measurements are stored in `test/validation/column_before_6197034.tsv` and
`column_after_e0ba262.tsv`. `exp/BettsMiller_Validation/columns.jl` produces
them from the fixed fixtures and analytic dry input; running that script with
the appropriate `--project` selects the implementation being measured.

## Five-day moist model comparison

`exp/BettsMiller_Validation/validate.jl` runs T21L20 on a 64 by 32 Gaussian
grid with even-sigma levels, from the balanced `:Moist_Spinup` initial state.
The pair uses identical initial fields and five simulated days at 600 s and
300 s. Both enable BM `:isca`, large-scale condensation with full heating,
surface sensible heat and evaporation, implicit PBL mixing, and Held--Suarez
forcing. The core's mass, dry-energy, and water corrections remain enabled;
the LRF is disabled. RH is 0.8 and tau is 7200 s.

The SST is explicitly prescribed as

```math
T_s(\phi)=271+29\exp\left[-\frac{\phi^2}{2(26\pi/180)^2}\right]\ \mathrm{K}.
```

This makes the comparison independent of the working-tree default SST
amplitude. The experiment writes hourly observations, final fields, and a
machine-readable result into the requested output directories. Its source
hashes identify the exact code used, including pre-existing working-tree edits.

### Budget method and scope

The production dynamics and coupled physics kernels advance the model. A
separate replay of the public physics processes observes their individual
increments and must match production states and averaged diagnostics at every
step. It records BM, condensation, sensible heat, evaporation, PBL mixing,
Rayleigh damping with heating, Newtonian relaxation, humidity cleanup, and
pressure adjustment. The normal spectral synchronization is measured too.
No convection formula is reimplemented in the observer.

The area-mean ledger tracks water, dry mass, and
$c_pT+L_vq+(u^2+v^2)/2$, using the current layer masses. Its cumulative history
follows the model's previous-level leapfrog recurrence; adding a `2Δt` flux
over every output sample would double-count the mature integration interval.
Each provisional dynamics state is compared with its previous-level source.
Dynamics, hyperdiffusion, time filtering, and its global fixers are grouped as
one measured net term. This grouping does not identify their individual
unfixed errors.

At fixed masses, BM and fully heated condensation conserve the moist-enthalpy
part of this ledger. Evaporation adds latent energy; sensible heat and
Newtonian relaxation add or remove thermal energy. Pressure adjustment can
change the pressure-weighted thermal/kinetic amount while preserving water
and column dry mass. PBL and pressure effects are measured separately, not
assumed to vanish. The ledger closes by accounting for these effects; it is
not a proof that this quantity is conserved by the complete model.

### Measured results, 2026-10-05

Both runs completed with Julia 1.8.5 on one thread. The initial fields match
bit for bit, and the source hashes match. Every observed replay matches
production physics exactly. Both deep and shallow BM adjustment occur, and
large-scale condensation and evaporation are nonzero. The initial area-mean
column water is `28.41441 kg/m²`.

| Quantity | 600 s | 300 s |
|:---|---:|---:|
| Completed steps / simulated days | 720 / 5 | 1440 / 5 |
| Temperature range over all steps, K | 168.474--305.200 | 168.460--305.200 |
| Specific humidity range, kg/kg | 0--0.0200790 | 0--0.0200978 |
| Accumulated BM rain, mm | 7.96371 | 7.89631 |
| Accumulated large-scale rain, mm | 1.32236 | 1.38551 |
| Accumulated evaporation, mm | 10.80302 | 10.82022 |
| Final area-mean column water, kg/m² | 29.93137 | 29.95281 |
| Largest BM temperature increment per substep, K | 0.5560 | 0.2792 |
| Largest BM temperature increment at an observed regime change, K | 0.4623 | 0.2367 |
| Largest BM humidity increment, kg/kg | 2.920e-4 | 1.421e-4 |

Halving the timestep changes accumulated BM rain by 0.85%, large-scale rain
by 4.56%, and evaporation by 0.16%, relative to the 300 s run. The two rain
components compensate: total accumulated rain differs by about 0.046%.
Final mass-weighted RMS differences are `0.2563 K` in temperature and
`1.4714e-4 kg/kg` in humidity. Local maxima are larger: `2.6977 K` and
`2.7743e-3 kg/kg`, respectively. Final wind RMS differences are `0.6690 m/s`
in u and `0.5457 m/s` in v; local differences reach about 9 m/s. This is
measured timestep sensitivity, not a claim of converged instantaneous fields.
Observed BM increments and increments at regime changes approximately halve
with the timestep, with no explosive jump or invalid state in this case.

The final water change agrees with accumulated evaporation minus both rain
components to `1.5e-12 kg/m²`. Maximum water-ledger residuals throughout the
runs are below `5.4e-12 kg/m²`; dry-air pressure adjustment residuals are below
`3.2e-10 kg/m²`. Maximum fixed-mass BM and condensation energy residuals,
measured by subtracting whole-state integrals, are below `4.8e-5 J/m²` per
substep. The complete energy ledger closes within `9.0e-4 J/m²` throughout
the run, relative to an initial ledger of about `2.62e9 J/m²`.

The measured five-day mean contributions to the energy ledger are:

| Process, W/m² | 600 s | 300 s |
|:---|---:|---:|
| Newtonian relaxation | -76.2173 | -76.2033 |
| Surface latent-energy input | 62.5175 | 62.6170 |
| Surface sensible heat | 6.1395 | 6.2922 |
| PBL mixing, measured net effect | 0.5834 | 0.5821 |
| Pressure adjustment, measured net effect | 0.8726 | 0.8853 |
| Net ledger change | -6.1043 | -5.8267 |

BM, fully heated condensation, and Rayleigh damping with heating have
negligible net contributions at this precision. The grouped corrected
dynamics and spectral-synchronization terms are also negligible in the global
ledger. PBL and pressure adjustment are **nonzero** for this ledger; omitting
them would produce an unexplained energy residual of about `1.46 W/m²`.
Their individual physical energy interpretation requires a broader model
study. No tuning or changes to those processes were made for this comparison.

The machine-readable comparison is
`test/validation/model_results_20261005.toml`; hourly observations are
`test/validation/model_hourly_{600,300}_20261005.tsv`. The report retains
per-process increments, extrema, replay errors, source hashes, and initial
state hashes. These results establish integration stability and explicit
budget accounting for the prescribed five-day case.

## Reproduction

From the repository root:

```sh
julia --project=. exp/BettsMiller_Validation/validate.jl /tmp/bm600 5 600
julia --project=. exp/BettsMiller_Validation/validate.jl /tmp/bm300 5 300
julia --project=. exp/BettsMiller_Validation/compare.jl \
    /tmp/bm600 /tmp/bm300 /tmp/bm-comparison.toml
```

Use the same Julia version and thread count for the paired comparison. The
comparison checks initial-state hashes and matching source hashes. Raw field
files belong in the chosen output folders; the small result report can be
retained with the validation record. The regular `exp/BettsMiller_Test/`
example also combines BM with condensation and writes the process diagnostics
through the normal output manager.

This is a prescribed short integration check. A climatology, validation of
each numerical-core correction, and a longer T21/T42 sensitivity study remain
separate scientific tasks. The corrections require rerunning previous model
experiments; choosing the old deep timescale option is not a reconstruction
of the old scheme.
