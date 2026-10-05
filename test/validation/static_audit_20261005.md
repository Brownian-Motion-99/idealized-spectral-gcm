# Static audit of Betts–Miller, 2026-10-05

Audited implementation: `e0ba262`, with the pre-existing working-tree changes
recorded separately before this audit. Reference: local Isca revision
`5c6d5baf74498c5dce3578e235c5254e6856baf6`,
`src/atmos_param/qe_moist_convection/qe_moist_convection.F90`. This is Frierson's
Simple Betts–Miller scheme with the O'Gorman–Schneider modifications described
in the module header.

The review covered `Betts_Miller.jl`, shared moist thermodynamics, driver
configuration, `Moist_Physics!`, `Spectral_Physics!`, dry-air adjustment, and
the dynamics-to-physics synchronization. No additional error in the ordinary
deep or shallow adjustment algebra was found. Step 5's independent Fortran
results support that conclusion; this audit also checks the control flow and
the interpretation of the returned diagnostics.

| Area | Static finding |
|:---|:---|
| Parcel ascent | Unsaturated ascent conserves actual parcel mixing ratio. Saturation mixing ratio is separate. Log-pressure LCL bisection brackets represented condensation and handles coincidence with full levels; moist integration only ascends. No represented condensation returns dry diagnostics and no adjustment. |
| Buoyancy | Parcel and environment use the dynamics' linear-in-specific-humidity virtual-temperature relation. The first contiguous moist buoyant region is retained; neutral roundoff does not initiate an LFC. Actual temperature enters ascent and enthalpy closure. |
| RH | Recovering vapor pressure from saturation mixing ratio, scaling it by RH, and converting to specific humidity gives the exact requested preliminary RH. It is intentionally not recomputed after the thermal closure. |
| Deep closure | `Pq=-sum(qdot*m)` and `PT=(cp/Lv)*sum(Tdot*m)` have consistent signs and precipitation units. When `Pq>PT`, only humidity rates are reduced. Otherwise the Isca uniform thermal correction cancels the moist-enthalpy residual. The alternative timescale mode reduces thermal rates. |
| Shallow closure | Removing upper layers uses surface-up moisture suffix sums. The accepted fractional cell satisfies `f*Qtop+Qlower=0`. Fractional rates are formed once, then full cell masses enter the uniform correction. An infeasible or single-cell adjustment is reset. Zero-rain transport remains active. |
| Boundaries | Exact zero integrals, interface crossings, out-of-domain LCLs, neutral crossings, and buoyancy reaching the model top have conservative explicit handling. These are deliberate improvements over source-reference edge behavior, not bitwise Isca equivalence. |
| Coupling | BM rates update the provisional next state once per bounded substep. Condensation sees that state. Diagnostics average substep rates. Pressure adjustment preserves column dry mass and each layer's post-physics water mass; subsequent synchronization restores global integrals after spectral projection. |
| Storage | Scalar calls own their buffers. Grid calls overwrite outputs and use distinct per-thread work arrays; scratch profiles and tendencies are reset for each column. A state should belong to one grid calculation at a time, with separate input/output arrays. |

Three interpretation and domain caveats need explicit documentation:

1. **CIN is an Isca-style diagnostic.** Below an in-domain LCL, the signed
   parcel/environment difference is integrated, so buoyant dry layers can
   reduce the diagnostic. Some inactive paths reset CIN to zero; the absent-LCL
   dry path instead accumulates only negative buoyancy. It is not universally
   a positive-only integral of all inhibition in the column. It does not
   enter the adjustment trigger.
2. **Standalone negative humidity is compatibility handling.** The column
   routine accepts finite `q<1` and evaluates humidity relative to `max(q,0)`
   without changing the caller. Its existing regression deliberately covers
   a negative spectral undershoot. Physical water budgets therefore apply to
   admissible nonnegative inputs, or to that cleaned input, not to arbitrary
   negative vapor mass. Coupled physics validates the resulting state and
   rejects material negative humidity. This compatibility behavior should not
   be described as accepting only roundoff-scale negatives.
3. **Column closure is not a whole-model energy theorem.** The BM identity is
   for `cp*T+Lv*q` at fixed pressure masses. Pressure changes, other physics,
   spectral projection, and numerical corrections must be accounted for in
   a model budget. `dt<=tau` bounds the explicit relaxation fraction, but does
   not independently guarantee positive temperatures or stability on every
   pressure grid. The moist RK2 discretization follows Isca and is not adaptive;
   unusually coarse vertical grids require additional accuracy checks.

Documentation and example follow-up: identify the reference formulation,
explain these diagnostic/domain limits, distinguish preliminary deep targets
from effective shallow targets, state that the previous results require
rerunning after the corrections, and remove the stale example claim that BM
and large-scale condensation cannot be combined. Step 6 should measure the
full model's sources, pressure-adjustment effects, and numerical residuals
before making an integration-stability claim.

The static audit is complete. No convection-physics change is required before
the step 6 experiments; the next work is documentation and measured model
validation.
