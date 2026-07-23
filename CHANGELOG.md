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

- The last timestep of a run no longer opens a new, empty trailing NC chunk.
  Previously, when `output_interval` divided `end_time` evenly (the normal
  case), the final flush wrote into the current chunk and then immediately
  rotated to a fresh chunk file that the run had no more data left to write,
  leaving a valid-but-empty `output_t<end_time>.nc` (correct header, zero
  time steps). `Driver.jl` now skips rotation once `integrator.time >=
  config.end_time`.
- `exp/LscaleCond_Test/plot_precip.py` updated for the chunked output format
  above. `netCDF4.MFDataset` cannot join the chunks — it only supports the
  `NETCDF3_*`/`NETCDF4_CLASSIC` formats, and `NCDatasets.jl` writes full
  `NETCDF4` (HDF5) — so the script now uses `xarray.open_mfdataset`, which
  also sorts chunks by their actual `time` coordinate rather than relying on
  filename order (important once time values span different digit counts,
  e.g. `output_t8640000.nc` vs `output_t17280000.nc`).

## [0.1.0]

Baseline prior to the chunked-output refactor.
