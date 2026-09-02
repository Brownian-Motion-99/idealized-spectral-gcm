# Model configuration

A simulation is configured with the keyword struct `Model_Config`. It contains
the equation set, grid, integration controls, initialization, I/O, and a
`Dict{String,Any}` of parameterized-physics options. `JGCM_Simulate` validates
cross-field constraints before allocating the model.

The complete type definition is in `src/Core/Experiment_Configuration.jl`.
Fields marked “required” below have no constructor default and must be passed
even when a particular equation set does not use them.

## Experiment identity and equation set

| Field | Type | Default | Meaning |
|:---|:---|:---|:---|
| `name` | `String` | required | experiment identifier used in logs and NetCDF metadata |
| `institution` | `String` | `"Unknown Institution"` | NetCDF institution metadata |
| `model_type` | `Symbol` | required | `:PrimitiveEquation`, `:Shallow_Water`, or `:Barotropic` |

The three values select different governing equations, not different physics
packages. Primitive-equation parameterizations are not called by the
barotropic or shallow-water drivers.

## Resolution and planetary constants

| Field | Type | Default | Meaning |
|:---|:---|:---|:---|
| `num_fourier` | `Int` | required | maximum configured zonal wavenumber $M$ |
| `nθ` | `Int` | required | number of Gaussian latitudes $n_\phi$ |
| `nd` | `Int` | required | number of vertical layers; use 1 for the 2-D models |
| `radius` | `Float64` | required | planetary radius $a$, m |
| `omega` | `Float64` | required | rotation rate $\Omega$, rad s$^{-1}$ |
| `grav` | `Float64` | required | gravitational acceleration $g$, m s$^{-2}$ |

The driver sets `nλ = 2nθ` and `num_spherical = num_fourier + 1`. It requires
`nθ` to be positive and even and checks

```math
n_\phi\geq\left\lceil\frac{3M+1}{2}\right\rceil.
```

Common transform grids are T21 with `nθ = 32` and T42 with `nθ = 64`.

## Vertical coordinate

These fields are used by `:PrimitiveEquation`. They are still required by the
configuration struct for 2-D modes, where `nothing` or placeholder strings are
accepted because no `Vert_Coordinate` is built.

| Field | Type | Supported values |
|:---|:---|:---|
| `vert_coord_option` | `Any` | `"even_sigma"`, `"uneven_sigma"`, `"hybrid"`, `"simmons_and_burridge"`, `"mcm"`, `"v197"` |
| `vert_difference_option` | `Any` | currently `"simmons_and_burridge"` |
| `vert_ref_level_option` | `Any` | `"second_centered_wts"` or `"second_centered"` |

`vert_ref_level_option` is the historical configuration name for the interface
interpolation used by vertical advection. The detailed definitions and the
fixed-table level restrictions are in [Vertical discretization](@ref).

## Time integration and numerical controls

| Field | Type | Default | Meaning |
|:---|:---|:---|:---|
| `Δt` | `Int64` | required | base timestep in seconds |
| `end_time` | `Int64` | required | duration advanced by this invocation, seconds |
| `day_to_sec` | `Int64` | required | seconds in the configured model day |
| `damping_order` | `Int` | required | even exponent of scale-selective damping |
| `damping_coef` | `Float64` | required | normalized grid-scale damping rate, s$^{-1}$ |
| `robert_coef` | `Float64` | required | Robert--Asselin filter coefficient |
| `implicit_coef` | `Float64` | required | semi-implicit weight $\alpha$ |

`end_time` must be positive and exactly divisible by `Δt`. On a warm start it
is still a **duration**, not an absolute target time: a checkpoint at day 100
with `end_time = 20day_to_sec` finishes at day 120.

The first integration call covers `Δt`; mature leapfrog calls connect time
levels separated by `2Δt`, while the model clock still advances by one base
timestep per driver iteration. The common centered setting is
`implicit_coef = 0.5`.

## Composition and initialization

| Field | Type | Default | Meaning |
|:---|:---|:---|:---|
| `moisture_processes` | `Bool` | `true` | transport and retain grid specific humidity |
| `initial_condition` | `Any` | required | built-in symbol, callback, or array mapping |

### Built-in initial conditions

| Symbol | Intended model | Description |
|:---|:---|:---|
| `:Moist_Spinup` | `:PrimitiveEquation` | Ullrich balanced shallow-atmosphere state with a localized wind perturbation and optional analytic moisture |
| `:Shallow_Water_Test` | `:Shallow_Water` | resting layer plus prescribed localized and ITCZ-like equilibrium-height anomalies |
| `:Barotropic_Jet` | `:Barotropic` | zonal jet with a localized wavenumber-4 vorticity perturbation |

`:Moist_Spinup` is also suitable for a dry run when
`moisture_processes = false`; its analytic virtual temperature then equals
actual temperature. When moisture is enabled, the optional physics dictionary
key `"initial_humidity_floor"` sets a nonnegative floor below 100 hPa.

### Custom callback

A function-valued `initial_condition` is called as

```julia
initial_condition(mesh, atmo_data, dyn_data, vert_coord)
```

It must populate a mutually consistent current state. This is a low-level
interface: a callback that sets grid winds or dry prognostics must also create
the corresponding spectral representation and initialize any required history.
Use the built-in initialization routines as templates.

### Arrays

A `Dict` or `NamedTuple` can inject arrays through the model-specific variable
map. Grid `:u` and `:v` are transformed together into spectral vorticity and
divergence. Grid `:t` and pressure/height aliases are transformed to their
spectral fields. Array shapes must exactly match allocated model fields. This
loader fills mapped current-level arrays only; it does not initialize every
previous-level history field used by primitive-equation conservation targets.
For a custom primitive-equation integration, prefer a callback that explicitly
initializes a consistent current and previous state.

## Numerical correction options

For historical reasons these dynamical-core controls are stored inside
`physics_params`:

| Key | Driver fallback | Meaning |
|:---|:---|:---|
| `"do_mass_correction"` | `true` | restore global mean surface pressure after dynamics |
| `"do_energy_correction"` | `true` | restore global pressure-mass-weighted $c_pT+K$ |
| `"do_water_correction"` | `true` | restore global grid-tracer water integral |
| `"use_virtual_temperature"` | `true` | use $T_v$ in moist equation-of-state terms |

These are numerical/core options, not parameterized physical processes. Moist
primitive-equation runs require `use_virtual_temperature = true`.

## Physical parameter dictionary

All primitive-equation schemes default to disabled when their master switch is
absent. Once a scheme is enabled, keys accessed directly by that scheme must be
provided.

| Process | Master switch | Other keys |
|:---|:---|:---|
| Held--Suarez | `"do_HS_Forcing"` | `"σ_b"`, `"k_a"`, `"k_s"`, `"k_f"`, `"ΔT_y"`, `"Δθ_z"`; optional `"T_equator"`, `"T_stratosphere"` |
| Betts--Miller | `"do_Betts_Miller"` | optional `"bm_tau"`, `"bm_relative_humidity"` |
| large-scale condensation | `"do_Lscale_Cond"` | optional `"condensation_heating_fraction"` or legacy `"L"` |
| sensible surface heat | `"do_Sensible_Heating"` | `"C_H"`; optional `"lower_boundary_temperature"` |
| surface evaporation | `"do_Surface_Evaporation"` | `"C_E"`; optional `"lower_boundary_temperature"` |
| PBL scalar mixing | `"do_Implicit_PBL_Scheme"` | `"C_D"`, `"PBL_Top_Mode"`, `"PBL_Top_Value"` |
| moisture LRF | `"do_LRF"` | `"LRF_file"` |

See [Physical parameterizations](@ref) for equations, units, validation, and
process order. A robust explicit dictionary includes every master switch, even
when false, so experiment intent is visible in the script.

## Shallow-water forcing keys

The shallow-water driver reads its idealized forcing and initialization values
from `physics_params`:

| Key | Fallback | Meaning |
|:---|:---|:---|
| `"h_0"` | `3.0e4` | reference layer geopotential used by the equations |
| `"kappa_m"` | `1/(20day_to_sec)` | linear momentum damping rate, s$^{-1}$ |
| `"kappa_t"` | `1/(10day_to_sec)` | equilibrium-height relaxation rate, s$^{-1}$ |
| `"h_amp"` | `2.0e4` | localized equilibrium-geopotential amplitude |
| `"h_lon"` | `90.0` | anomaly longitude, degrees |
| `"h_lat"` | `25.0` | anomaly latitude, degrees |
| `"h_width"` | `15.0` | anomaly width, degrees |
| `"h_itcz"` | `1.0e5` | equatorial equilibrium-geopotential amplitude |
| `"itcz_width"` | `4.0` | equatorial feature width, degrees |

`"alpha"` appears in the supplied shallow-water experiment dictionary but is
not read by the current implementation.

## Restart controls

| Field | Type | Default | Meaning |
|:---|:---|:---|:---|
| `is_restart` | `Bool` | `false` | load a checkpoint instead of applying `initial_condition` |
| `restart_file` | `String` | `""` | JLD2 checkpoint path for a warm start |
| `saving_frequency` | `Int64` | `86400` | checkpoint and NetCDF chunk interval in seconds; 0 disables both |

If positive, `saving_frequency` must be divisible by `Δt`. Checkpoint behavior,
retention, and cold-start cleanup are detailed in [Output and restarts](@ref).

## Output controls

| Field | Type | Default | Meaning |
|:---|:---|:---|:---|
| `output_path` | `String` | required | base directory and restart-directory parent |
| `output_filename` | `String` | required | base path used to derive NetCDF chunk names |
| `logger` | `String` | required | runtime log path |
| `vars_to_output` | `Vector{Symbol}` | required | requested model-specific diagnostics |
| `output_interval` | `Int64` | required | averaging and flush interval, seconds |
| `spinup_day` | `Float64` | `0.0` | time after segment start during which samples are discarded, model days |
| `do_plev_output` | `Bool` | `false` | additionally write pressure-level output |
| `pressure_levels` | `Vector{Float64}` | standard ten-level vector | requested pressure levels in Pa |

`output_interval` must be positive and divisible by `Δt`. Pressure-level output
is available only for `:PrimitiveEquation`, requires a nonempty vector of
positive finite pressures, and always adds `:ps` to the active variable list.
See [Output and restarts](@ref) for the variable namelist and file semantics.

## Validation summary

Before allocation the driver rejects:

- unknown `model_type` values;
- negative truncation, nonpositive or odd `nθ`, an under-resolved transform
  grid, or nonpositive `nd`;
- nonpositive `radius`, `grav`, `day_to_sec`, `Δt`, `end_time`, or
  `output_interval`;
- durations or output/checkpoint intervals not divisible by `Δt`;
- negative `saving_frequency` or `spinup_day`;
- pressure-level output for a 2-D model or with invalid target pressures;
- moist primitive-equation runs without virtual temperature;
- Betts--Miller runs with `Δt > bm_tau`;
- Betts--Miller, condensation, or LRF without moisture transport.

Individual parameterizations perform additional checks when called, including
explicit Held--Suarez stability bounds, humidity/temperature validity, PBL top
types, saturation-domain checks, and external LRF array dimensions.
