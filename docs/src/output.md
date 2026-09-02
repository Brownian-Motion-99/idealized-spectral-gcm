# Output and restarts

The driver writes time-mean NetCDF diagnostics, a plain-text runtime log, and
optional JLD2 checkpoints. NetCDF output and checkpoints share the same model
clock but have independent intervals: `output_interval` controls diagnostic
records and `saving_frequency` controls checkpointing and file rotation.

## NetCDF file naming

`output_filename` is a base name. If it is `run/output.nc` and the segment
starts at time zero, the native-grid file is

```text
run/output_t0.nc
```

When pressure-level output is enabled, the paired file is

```text
run/output_t0_plev.nc
```

At each positive `saving_frequency` boundary the driver writes a checkpoint,
closes the current NetCDF file, and opens a new chunk whose suffix is the
absolute model time in seconds. It does not open an empty new chunk at the
final boundary of a segment.

## Temporal averaging

After every base timestep, active variables are accumulated. When absolute
time relative to the current segment start is a multiple of
`output_interval`, the arithmetic mean of the accumulated samples is written.
The time coordinate is the record's ending model time in units of model days.

The initial state at the exact segment start is not accumulated. With
`spinup_day > 0`, all samples through

```math
t_{spinup}=t_{segment\ start}+
\mathtt{spinup\_day}\times\mathtt{day\_to\_sec}
```

are skipped. Finalization and chunk rotation flush a partially accumulated
interval, so its record can contain fewer samples than a regular interval.

For online pressure output, native-coordinate fields are averaged first and
the time-mean field is then interpolated. This is not the same as averaging
instantaneously interpolated fields.

## Primitive-equation variables

The following symbols are accepted in `vars_to_output`. `:ps` is added
automatically because native-coordinate and pressure-level metadata require
surface pressure.

### State and diagnostic fields

| Symbol | NetCDF name | Long name | Units | Dimensions |
|:---|:---|:---|:---|:---|
| `:u` | `ua` | eastward wind | m s$^{-1}$ | 3-D |
| `:v` | `va` | northward wind | m s$^{-1}$ | 3-D |
| `:w` | `wap` | pressure vertical velocity $\omega$ | Pa s$^{-1}$ | 3-D |
| `:t` | `ta` | air temperature | K | 3-D |
| `:q` | `hus` | specific humidity | 1 | 3-D |
| `:vor` | `vor` | relative vorticity $\zeta$ | s$^{-1}$ | 3-D |
| `:div` | `div` | divergence $\delta$ | s$^{-1}$ | 3-D |
| `:p` | `p` | full-level air pressure | Pa | 3-D |
| `:z` | `zg` | geopotential height | m | 3-D |
| `:ps` | `ps` | surface pressure | Pa | 2-D |
| `:lnps` | `lnps` | log surface pressure | 1 | 2-D |

Requesting `:z` also requires `:t`, because pressure-level geopotential-height
interpolation uses temperature hydrostatically.

### Physical-process diagnostics

| Symbol | NetCDF name | Long name | Units | Dimensions |
|:---|:---|:---|:---|:---|
| `:t_eq` | `teq` | Held--Suarez equilibrium temperature | K | 3-D |
| `:shflx` | `hfss` | upward surface sensible heat flux | W m$^{-2}$ | 2-D |
| `:lhflx` | `hfls` | upward surface latent heat flux | W m$^{-2}$ | 2-D |
| `:precip` | `pr` | total precipitation flux | kg m$^{-2}$ s$^{-1}$ | 2-D |
| `:bm_dt` | `bm_dta_dt` | Betts--Miller temperature tendency | K s$^{-1}$ | 3-D |
| `:bm_dq` | `bm_dhus_dt` | Betts--Miller humidity tendency | s$^{-1}$ | 3-D |
| `:bm_precip` | `bm_pr` | Betts--Miller precipitation flux | kg m$^{-2}$ s$^{-1}$ | 2-D |
| `:lrf_dt` | `lrf_dta_dt` | LRF temperature tendency | K s$^{-1}$ | 3-D |

### Model tendency fields

| Symbol | NetCDF name | Long name | Units | Dimensions |
|:---|:---|:---|:---|:---|
| `:du` | `dua_dt` | eastward-wind tendency | m s$^{-2}$ | 3-D |
| `:dv` | `dva_dt` | northward-wind tendency | m s$^{-2}$ | 3-D |
| `:dt` | `dta_dt` | temperature tendency | K s$^{-1}$ | 3-D |
| `:dq` | `dq_dt` | specific-humidity tendency | s$^{-1}$ | 3-D |
| `:dvor` | `dvor_dt` | relative-vorticity tendency | s$^{-2}$ | 3-D |
| `:ddiv` | `ddiv_dt` | divergence tendency | s$^{-2}$ | 3-D |
| `:dps` | `dps_dt` | surface-pressure tendency | Pa s$^{-1}$ | 2-D |

These are the model's working tendency diagnostics at output time; they should
not be assumed to be a complete process-budget decomposition. Use the named
Betts--Miller and LRF fields when those isolated process rates are required.

## Barotropic variables

| Symbol | NetCDF name | Units |
|:---|:---|:---|
| `:vor` | `vor` | s$^{-1}$ |
| `:u` | `ua` | m s$^{-1}$ |
| `:v` | `va` | m s$^{-1}$ |
| `:ke` | `ke` | m$^2$ s$^{-2}$ |
| `:dvor` | `dvor` | s$^{-2}$ |

All barotropic variables are two-dimensional.

## Shallow-water variables

| Symbol | NetCDF name | Units |
|:---|:---|:---|
| `:h` | `height` | m in current metadata; see note below |
| `:u` | `ua` | m s$^{-1}$ |
| `:v` | `va` | m s$^{-1}$ |
| `:vor` | `vor` | s$^{-1}$ |
| `:div` | `div` | s$^{-1}$ |
| `:pv` | `pv` | m$^{-1}$ s$^{-1}$ |
| `:dh` | `dh_dt` | m s$^{-1}$ |

All shallow-water variables are two-dimensional.

The shallow-water momentum equation adds `:h` to specific kinetic energy, so
the numerical field behaves as geopotential (m$^2$ s$^{-2}$), despite the
current output metadata calling it geopotential height in metres. The metadata
and the documented `:pv` units inherit that inconsistency.

## Native vertical-coordinate files

Primitive-equation native files use dimensions `(lon, lat, pfull, time)` for
three-dimensional fields. They include Isca-style vertical metadata:

- `pfull` and `phalf` model-level coordinates;
- `pk` and `bk`, the interface $A$ and $B$ coefficients;
- `ps`, surface pressure.

Full-level pressure is computed with the Simmons--Burridge layer-mean log
pressure, not with arithmetic midpoint coefficients. Generic CF readers cannot
fully reconstruct that nonlinear definition from the current metadata alone.
Use the model's pressure-level output or supplied postprocessor when exact
vertical placement matters.

## Pressure-level interpolation

Set

```julia
do_plev_output = true
pressure_levels = [92500.0, 85000.0, 50000.0, 20000.0]
```

to create a paired pressure-level file. Target values are in Pa and may be in
any requested order. Pressure output is available only for
`:PrimitiveEquation`.

For each time-mean record, the interpolation path:

1. reconstructs full-level log pressure from time-mean $p_s$ and the exact
   interface coefficients;
2. searches for the two bracketing full levels once per target point;
3. reuses those indices and log-pressure weights for every 3-D field;
4. uses a hydrostatic formula for `:z` and linear interpolation in log pressure
   for other fields.

Targets below the lowest model full level are written as `NaN`. Targets above
the top full level use a limited extrapolation weight clamped to
`[-0.5, 1.5]`. Two-dimensional fields such as $p_s$ are copied without vertical
interpolation.

### Offline interpolation

For long integrations, writing only native output is usually faster. Convert a
completed native file with:

```shell
julia --project=. post_processing/Interpolator.jl \
    input.nc output_plev.nc 92500 85000 50000 20000
```

The postprocessor reads `pk`, `bk`, and `ps`, uses the same cached interpolation
implementation as online output, copies the existing time coordinate, and
preserves variable attributes. Its default 3-D variable set is `:u`, `:v`,
`:t`, `:z`, `:q`, `:vor`, and `:div`; fields absent from the input are skipped.
The output path is replaced if it already exists.

## Checkpoints and warm starts

With `saving_frequency > 0`, checkpoints are written atomically as

```text
<output_path>/restart/restart_t<seconds>.jld2
```

A temporary file is completed and then renamed into place. Versioned
checkpoints store the absolute saved time, resolution tuple, and every array in
`Dyn_Data`, including all three time levels required to resume leapfrog
integration. Resolution is checked on load.

### Cold start

When `is_restart = false`, `initial_condition` is applied and model time starts
at zero. If `<output_path>/restart/` already exists, the driver deletes every
`.jld2` file in that directory before initialization. Use a new output
directory or copy checkpoints elsewhere before launching a cold start when old
files must be retained.

### Warm start

Set:

```julia
is_restart = true
restart_file = "run/restart/restart_t86400.jld2"
```

The saved state and absolute time are restored, the Euler startup is skipped,
and the semi-implicit matrices are rebuilt for the mature leapfrog interval.
`end_time` specifies how much **additional** time to integrate.

At checkpoint boundaries, the driver retains the latest five restart files.
During a warm-start segment it protects the checkpoint used to start that
segment from cleanup even if it is older.

## Runtime log

`logger` receives initialization information, checkpoint messages, periodic
state diagnostics, and a final performance summary. Progress is reported at
approximately quarter-model-day intervals and includes:

- segment step and percentage;
- absolute model day;
- elapsed wall time and ETA for the current invocation;
- maxima and locations for selected state and diagnostic fields;
- final seconds per step and simulated model days per wall day.

The log is overwritten at the start of each invocation, so preserve it before
launching another segment if the prior run history is needed.
