# Changelog

All notable changes to this project are documented in this file.

## [0.3.0]

### Changed (breaking)

- **`Model_Config` field renamed:** `restart_frequency` → `saving_frequency`. This
  single value now controls both JLD2 checkpoint cadence and NC chunk rotation —
  they always happen together. Scripts still passing `restart_frequency` will
  fail immediately with a `MethodError`.
- **`Model_Config` field renamed and semantics changed:** `do_raw_output` →
  `do_plev_output`. Native sigma/hybrid-sigma output is now *always* written
  (previously it was the `do_raw_output=true` case); `do_plev_output` only
  controls whether an *additional* pressure-level-interpolated file is also
  produced. There is no longer a way to get pressure-level output only.
- **Output files are always time-stamped chunks**, not a single monolithic
  file: `output_t0.nc`, `output_t<N>.nc`, ... (and `output_t<N>_plev.nc` when
  `do_plev_output=true`), rotating at each `saving_frequency` boundary. Any
  downstream script assuming a single `output.nc` needs to glob for
  `output_t*.nc` instead. (`exp/LscaleCond_Test/plot_precip.py` updated
  accordingly, using `xarray.open_mfdataset` to join chunks — see note below.)

### Fixed

- **Warm-restart chunk rotation** (`Driver.jl`): The `Rotate_NC_Chunk!` skip
  condition was comparing `integrator.time` against `config.end_time` (the
  *duration* of the invocation). On a warm restart, `integrator.time` starts
  at the checkpoint time and is therefore always greater than `config.end_time`
  for any non-trivial prior run, silently disabling all mid-run chunk rotation.
  The condition is now `integrator.time >= segment_end_time`
  (`= start_time + config.end_time`, the absolute model time at the end of
  the current invocation), which works correctly for both cold and warm starts.
- **No empty trailing NC chunk**: `Driver.jl` skips rotation on the final
  saving boundary of the run (when `integrator.time >= segment_end_time`),
  preventing a valid-but-empty `output_t<segment_end_time>.nc` file from being
  created after the last data has already been flushed into the previous chunk.
- `exp/LscaleCond_Test/plot_precip.py` updated for the chunked output format.
  `netCDF4.MFDataset` cannot join the chunks — it only supports the
  `NETCDF3_*`/`NETCDF4_CLASSIC` formats, and `NCDatasets.jl` writes full
  `NETCDF4` (HDF5) — so the script now uses `xarray.open_mfdataset`, which
  also sorts chunks by their actual `time` coordinate rather than relying on
  filename order (important once time values span different digit counts,
  e.g. `output_t8640000.nc` vs `output_t17280000.nc`).

### Tests

- Added `test/test_output_system.jl` Test 4: two-phase warm-start integration
  test. Phase 1 runs a 2-day cold start with `saving_frequency = 1 day` to
  produce a day-2 checkpoint. Phase 2 restarts from that checkpoint and runs
  for 1 day with `saving_frequency = 2 days` (no boundary hit). Asserts a
  single NC chunk at `t=172800`, correct 24-point time coverage (day 2+1h →
  day 3), no spurious trailing chunk, and no new JLD2 checkpoint.

## [0.1.0]

Baseline prior to the chunked-output refactor.
