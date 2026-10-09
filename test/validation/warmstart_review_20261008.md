# Warm-start validation, 2026-10-08

The configured T42L20 checkpoint resumes correctly with the current working
tree. All 74 persisted arrays load byte for byte with their original shapes
and element types. The accepted warm-start state and all four resumed steps
also match the independent continuation byte for byte. All 7,476 assertions
pass under Julia 1.8.5 with one thread.

The tested input is
`/data92/garywu/undergrad_proposal/ctrl_BM/restart/restart_t315360000.jld2`,
format version 3, with dimensions `(42, 43, 128, 64, 20)` and absolute saved
time `315360000` seconds. Its SHA-256 before and after testing is
`07cb1e4143f345be4fc77948d4ab05e810171babf7bada0bc2276c26b6fb522c`.
Generated checkpoints, NetCDF output, and runtime logs are in
`/tmp/jgcm_warmstart_review_20261008_final/`.

The field inventory covers every top-level array in `Dyn_Data`, including
all previous/current/next prognostic buffers, spectral tendencies, pressure
and geopotential diagnostics, surface fluxes, precipitation, BM and LRF
diagnostics, and scratch arrays:

| Field category | Shape | Count |
|:---|:---|---:|
| Spectral, full levels | `(43, 44, 20)` | 16 |
| Spectral, surface | `(43, 44, 1)` | 4 |
| Grid, full levels | `(128, 64, 20)` | 35 |
| Grid, surface | `(128, 64, 1)` | 13 |
| Grid, half levels | `(128, 64, 21)` | 6 |
| Total | | 74 |

Each saved value was read independently from JLD2 and compared against a
destination whose arrays had first been filled with NaNs. Every destination
buffer retained its identity. No saved fields were missing or nonfinite.
The loader validates the resolution tuple and rejects incorrect field shapes
on each of the three grid axes. Synthetic serialization tests fill every
element of every array with distinct values to detect omitted fields or
swapped time-level buffers. Existing version-2 and legacy compatibility tests
also pass.

Three driver comparisons run six steps continuously and compare against
three cold-start steps followed by three warm-start steps. The production
driver's state was observed immediately before its first resumed step, then
each of the three resumed checkpoints was compared against the continuous
run. Both experiments enable Held–Suarez forcing, condensation, surface
heating/evaporation, PBL mixing, and conservation corrections; the BM switch
varies as shown. LRF is disabled.

| Experiment | Loaded arrays | Arrays after each resumed step | Maximum absolute difference |
|:---|:---|:---|---:|
| T5L8, BM off, continuous versus split | 74/74 exact | 74/74 exact, 3 steps | 0 |
| T5L8, BM on, continuous versus split | 74/74 exact | 74/74 exact, 3 steps | 0 |
| T42L20, BM on, continuous versus split | 74/74 exact | 74/74 exact, 3 steps | 0 |
| Actual T42L20 checkpoint, BM on | 74/74 exact | 74/74 exact, 4 steps | 0 |

For the actual checkpoint, an independent in-memory continuation provides the
step-by-step reference. The test observes the real driver before dynamics
through a test-process method that delegates to its original dispatcher;
production code is unchanged. The semi-implicit wave matrix is checked
against the independently prepared leapfrog solver.

The clock checks verify the following at every driver step:

- Absolute `time` and `start_time` restore to `315360000` seconds.
- The configured timestep is `600` seconds, with a mature leapfrog interval
  of `1200` seconds. The Euler startup is skipped on every warm-start step.
- The segment duration is `2400` additional seconds, so `end_time` becomes
  `315362400` seconds.
- Checkpoints occur at `315360600`, `315361200`, `315361800`, and `315362400`
  seconds. NetCDF record times agree with these absolute interval-end times.
- A deliberately failing cold-initialization callback is never invoked.

Two implementation details bound the result. `Δt` comes from `Model_Config`;
the restart file stores absolute time but does not serialize or validate the
original timestep. The loader also skips absent array datasets for backward
compatibility, so it does not itself reject every incomplete checkpoint.
This audit independently checks completeness; the tested checkpoint contains
all 74 fields. Model identity and resolution remain configuration metadata,
and the nested tracer workspace is recreated as scratch storage.

Primitive-equation initialization reconstructs `grid_vor` and `grid_div` from
the restored current spectral fields to support older checkpoints with stale
diagnostics. In the final validation these two derived arrays differ from
their saved versions by at most `5.42e-20` and `1.36e-20`, respectively;
the other 72 arrays remain byte for byte unchanged. Both driver and reference
produce identical reconstructed diagnostics and identical resumed states.
The existing post-physics restart regression also passes its checks for a
checkpoint containing deliberately stale grid diagnostics.

Evidence is retained in
[the result and clock traces](warmstart_review_20261008_results.toml).
The field inventory and test results are summarized above; duplicate console
logs and the generated field table are omitted from the repository.
The result file includes SHA-256 hashes for every Julia source file and the
supporting test scripts. The strengthened all-field serialization checks are
in [test_Memory_And_Restart.jl](../test_Memory_And_Restart.jl).

To reproduce from the repository root with a new output directory:

```sh
julia --project=. --startup-file=no --threads=1 \
  test/validation/warmstart_review_20261008.jl \
  /tmp/jgcm-warmstart-new \
  /data92/garywu/undergrad_proposal/ctrl_BM/restart/restart_t315360000.jld2
```

Use `none` as the final argument to run only the generated-checkpoint
comparisons. The harness refuses to reuse an existing output directory.
