# Prescribed BM model validation

Run the five-day T21L20 moist comparison from the repository root:

```sh
julia --project=. exp/BettsMiller_Validation/validate.jl /tmp/bm600 5 600
julia --project=. exp/BettsMiller_Validation/validate.jl /tmp/bm300 5 300
julia --project=. exp/BettsMiller_Validation/compare.jl \
    /tmp/bm600 /tmp/bm300 /tmp/bm-comparison.toml
```

The setup enables Isca BM, fully heated large-scale condensation, surface
exchange, PBL mixing, and Held–Suarez forcing, with the normal core corrections.
An explicit zonally symmetric SST avoids dependence on edits to the default
surface temperature. Both runs use the same analytic initial fields; their
hashes and source hashes must match for the paired comparison.

`validate.jl` advances the production dynamics and coupled physics kernels.
A separate replay through the public physics processes records individual
budget changes and must match production states and averaged diagnostics at
every step. A previous-level leapfrog recurrence accounts for the budget and
precipitation histories. The energy ledger is `cp*T+Lv*q+K`, with PBL, pressure
adjustment, dynamics/fixers, and spectral synchronization explicitly measured.
Its closure is not a full-model energy-conservation theorem.

Each output folder contains `initial.bin`, `final.bin`, `hourly.tsv`, and
`result.toml`. Binary fields use Julia Serialization. Results for 2026-10-05
are retained in `test/validation/model_results_20261005.toml` and the two
`model_hourly_*_20261005.tsv` files; generated fields remain outside git.

See [the validation documentation](../../docs/src/betts-miller-validation.md)
for measured timestep sensitivity, physical interpretation, and limitations,
and [the prior static audit](../../test/validation/static_audit_20261005.md).
`columns.jl` additionally measures the fixed deep/shallow/virtual columns and
an analytic dry-neutral input; using a pre-fix project with the same script
reproduces the documented before/after diagnostics.
