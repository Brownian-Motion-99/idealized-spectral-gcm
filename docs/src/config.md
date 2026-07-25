# Model Configuration & Initialization

The `idealized-spectral-gcm` uses a two-tiered configuration system to separate **numerical/computational settings** (managed via a rigid `struct`) from **physical parameterizations** (managed via a flexible `Dict`).

## Physics Configuration (`physics_params`)

The `physics_params` dictionary toggles physical processes and sets their coefficients. These parameters are passed into the `Atmo_Data` structure during initialization or accessed directly by the forcing modules.

### Global Conservation Corrections

These flags enable global multiplicative correctors to enforce conservation laws, compensating for numerical drift or filter dissipation.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `"do_mass_correction"` | `Bool` | `true` | Rescales surface pressure to conserve total dry air mass. |
| `"do_energy_correction"` | `Bool` | `true` | Rescales global temperature to conserve total energy. |
| `"do_water_correction"` | `Bool` | `true` | Rescales specific humidity to conserve total moisture. |
| `"use_virtual_temperature"` | `Bool` | `true` | Uses virtual temperature ($T_v$) in hydrostatic balance and density calculations. |

### Large-scale Condensation

This is the configuration of `Lscale_Cond.jl`, controlling grid scale condensation.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `"do_Lscale_Cond"` | `Bool` | `true` | Master switch for the module. |
| `"L"` | `Float64` | `0.2` | Latent heating efficiency (ranges from 0.0 to 1.0). |

### Planetary Boundary Layer Processes

Configuration of `PBL.jl`, handling PBL-related processes, including surface sensible heat fluxes, surface latent heat fluxes, and PBL mixing for sensible/latent heat.

| Key | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `"do_Sensible_Heating"` | `Bool` | `true` | Toggle of surface sensible heat fluxes. |
| `"C_H"` | `Float64` | `0.0044` | Sensible heat flux coefficient. |
| `"do_Surface_Evaporation"` | `Bool` | `true` | Toggle of surface evaporation (latent heat fluxes). |
| `"C_E"` | `Float64` | `0.0044` | Latent heat flux (evaporation) coefficient. |
| `"do_Implicit_PBL_Scheme"` | `Bool` | `true` | Toggle of PBL mixing. |
| `"C_D"` | `Float64` | `0.0044` | Drag coefficient. |
| `"PBL_Top_Mode"` | `Symbol` | `:PressureLevel` | PBL top mode (`:PressureLevel` or `:ModelLevel`). |
| `"PBL_Top_Value"` | `Any` | `85000.0` | PBL top value (`Float64` or `int`, `Float64` for `:PressureLevel`, `int` for `:ModelLevel`). |


### Held-Suarez Forcing Parameters

Controlled by `HS_Forcing.jl`. These parameters define the relaxation toward the zonally symmetric equilibrium state.

| Key | Type | Default | Physical Meaning |
| :--- | :--- | :--- | :--- |
| `"do_HS_Forcing"` | `Bool` | `true` | Master switch for the module. |
| `"σ_b"` | `Float64` | `0.7` | Top of the planetary boundary layer (sigma coordinate). |
| `"k_f"` | `Float64` | `1.0` day$^{-1}$ | Surface friction damping rate. |
| `"k_a"` | `Float64` | `1/40` day$^{-1}$ | Thermal relaxation rate in the free atmosphere. |
| `"k_s"` | `Float64` | `1/4` day$^{-1}$ | Thermal relaxation rate at the surface. |
| `"ΔT_y"` | `Float64` | `60.0` K | Equator-to-pole temperature difference. |
| `"Δθ_z"` | `Float64` | `10.0` K | Vertical potential temperature gradient. |

> **Note on Energy Conservation:** The Held-Suarez module in this dycore explicitly calculates **frictional heating**. The kinetic energy dissipated by the Rayleigh damping term (`"k_f"`) is returned to the thermodynamic equation as sensible heat, ensuring closed energy budgets.

---

## Core Experiment Configuration (`Model_Config`)

The `Model_Config` struct (defined in `Experiment_Configuration.jl`) controls the grid generation, memory allocation, and time integration.

### Resolution & Geometry

| Parameter | Type | Description |
| :--- | :--- | :--- |
| `num_fourier` | `Int` | Spectral truncation limit ($M$). For T21, set to 21. |
| `nθ` | `Int` | Number of Gaussian latitudes ($J$). **Constraint:** $J \geq (3M+1)/2$. |
| `nd` | `Int` | Number of vertical levels. |
| `radius` | `Float64` | Planetary radius ($a$). Default: $6.371 \times 10^6$ m. |
| `omega` | `Float64` | Angular velocity ($\Omega$). Default: $7.292 \times 10^{-5}$ rad/s. |

### Vertical Coordinate System

The `vert_coord_option` parameter selects the strategy for generating interface levels ().

* **`"simmons_and_burridge"`**: The standard 20-level hybrid  coordinate used in NCAR CAM and the original spectral core.
* **`"even_sigma"`**: Uniformly spaced sigma levels from 0 to 1.
* **`"uneven_sigma"`**: Stretched sigma levels with higher resolution near the surface.
* *Requires:* `scale_heights` (stretching factor), `surf_res` (surface resolution), `exponent`.


* **`"hybrid"`**: A custom blend using a sine-squared transition function.
* *Requires:* `p_sigma` (pure sigma threshold), `p_press` (pure pressure threshold).



### Time Integration & Numerics

| Parameter | Standard Value | Description |
| :--- | :--- | :--- |
| `Δt` | `600` | Time step in seconds. |
| `end_time` | - | Duration of the current simulation invocation in seconds. |
| `damping_order` | `4` | Order of hyper-diffusion ($\nabla^4$). |
| `damping_coef` | `1.15e-4` | Diffusion coefficient. |
| `robert_coef` | `0.04` | Robert-Asselin time filter strength. |
| `implicit_coef` | `0.5` | Semi-implicit weighting (0.5 = Centered/Crank-Nicolson). |

---

## Initialization & Restart Strategy

The `Driver.jl` module handles the simulation lifecycle. The behavior depends heavily on the `is_restart` flag.

### Cold Start (`is_restart = false`)

Use this for new experiments.

* **Initialization:** Sets analytical initial conditions defined by `initial_condition` (e.g., `:Moist_Spinup`).

**WARNING:** The driver will **delete all existing `.jld2` files** in the `restart/` subdirectory to prevent mixing data from previous runs. Save the restart files in different directory is strongly recommended.

### Warm Start (`is_restart = true`)

Use this to continue a simulation.

* **Initialization:** Loads the full model state (spectral coefficients + surface fields) from `restart_file`.
* **Time:** The model clock resumes exactly where the restart file left off and advances for `end_time` seconds. For example, a restart at model day 10,000 with `end_time` set to 20,000 days finishes at model day 30,000.
* **Filenames:** New restart files and NC output chunks will be generated periodically based on `saving_frequency` (in seconds).

---

## Diagnostics & Output

### Standard Output (NetCDF)

The `Output_Manager` writes time-averaged snapshots to NetCDF files at `output_interval` (in seconds).

* **Fields:** Defined by `vars_to_output`. See [Available Output Namelist](#available-output-namelist) for available outputs.
* **Primary output:** Data is always written on the native sigma/hybrid-sigma coordinate grid (e.g., `output_t0.nc`).
* **Chunked files:** When `saving_frequency > 0`, the output is split into time-stamped chunks (`output_t0.nc`, `output_t86400.nc`, …) coordinated with JLD2 checkpoint saves. Each completed chunk is safe even if the run crashes later.
* **Pressure-level output (optional):** Set `do_plev_output = true` together with `pressure_levels` to also produce interpolated output on pressure levels (e.g., `output_t0_plev.nc`). An error is raised if `do_plev_output = true` but `pressure_levels` is empty.
* **Final chunk:** `Driver.jl` does not rotate to a new NC chunk on the final saving boundary of the run. The rotation is skipped when `integrator.time >= segment_end_time` (the absolute model time at the end of the current invocation, i.e. `start_time + end_time`). This correctly handles both cold starts (where `segment_end_time == end_time`) and warm restarts (where `segment_end_time > end_time`). The last interval's data stays in the previously-opened chunk rather than being moved to a fresh, empty file.

### Runtime Logging

The driver prints a status summary to `logger.log` periodically.

* **Monitoring:** It tracks the maximum absolute value and location of selected state and diagnostic fields, including $U$, $V$, $T$, vertical velocity, pressure, equilibrium temperature, and moisture.
* **Segment step:** A segment is the work performed by the current invocation. Both cold and warm starts advance for the duration specified by `end_time`; a warm start begins at the checkpoint's model time. For example, `Segment Step 360/10080` means that the current invocation has completed 360 of its 10,080 timesteps.
* **Progress:** `Segment` is the percentage completed by the current invocation. Elapsed time and ETA also apply to the current segment. ETA uses the average wall-clock time per completed segment step and becomes more accurate as the run progresses. The displayed model day remains the absolute model time, including time accumulated before a warm start.
* **Usage:** Use this log to follow run progress and detect numerical instability, such as rapidly growing field magnitudes.

At successful completion, the driver logs a `Run Metrics` summary containing:

* The model-day interval advanced by the current invocation.
* The number of timesteps completed.
* Initialization, integration, and total wall-clock time.
* Average wall-clock seconds per timestep.
* **Simulated days/wall day:** A throughput rate extrapolated to 24 hours of real runtime. For example, `5151.07 simulated days/wall day` means that, at the measured speed, one wall-clock day would advance the model by approximately 5,151 simulated days.

---

Based on the `Output_Mappings.jl` file you provided, here is the formatted **Available Output Namelist** section. You can append this to the end of your `configuration.md` or place it within the **Diagnostics & Output** section.

---

### Available Output Namelist

The following tables list the valid symbols available for the `vars_to_output` list in `Model_Config`. The availability depends on the `model_type` selected.

#### Primitive Equation (`:PrimitiveEquation`)

These variables cover full 3D atmospheric states and 2D surface fluxes.

| Symbol | NetCDF Name | Units | Standard Name (`std_name`) | Dimensions |
| :--- | :--- | :--- | :--- | :--- |
| **State Variables** | | | | |
| `:u` | `ua` | m s-1 | `eastward_wind` | 3D |
| `:v` | `va` | m s-1 | `northward_wind` | 3D |
| `:w` | `wap` | Pa s-1 | `lagrangian_tendency_of_air_pressure` | 3D |
| `:t` | `ta` | K | `air_temperature` | 3D |
| `:q` | `hus` | 1 | `specific_humidity` | 3D |
| `:vor` | `vor` | s-1 | `atmosphere_relative_vorticity` | 3D |
| `:div` | `div` | s-1 | `divergence_of_wind` | 3D |
| `:p` | `p` | Pa | `air_pressure` | 3D |
| `:z` | `zg` | m2 s-2 | `geopotential` | 3D |
| `:t_eq` | `teq` | K | `held_suarez_equilibrium_temperature` | 3D |
| **Surface & Fluxes** | | | | |
| `:ps` | `ps` | Pa | `surface_air_pressure` | 2D |
| `:lnps` | `lnps` | 1 | `log_surface_air_pressure` | 2D |
| `:shflx` | `hfss` | W m-2 | `surface_upward_sensible_heat_flux` | 2D |
| `:lhflx` | `hfls` | W m-2 | `surface_upward_latent_heat_flux` | 2D |
| `:precip` | `pr` | kg m-2 s-1 | `precipitation_flux` | 2D |
| **Tendencies** | | | | |
| `:du` | `dua_dt` | m s-2 | `tendency_of_eastward_wind` | 3D |
| `:dv` | `dva_dt` | m s-2 | `tendency_of_northward_wind` | 3D |
| `:dt` | `dta_dt` | K s-1 | `tendency_of_air_temperature` | 3D |
| `:dps` | `dps_dt` | Pa s-1 | `tendency_of_surface_air_pressure` | 2D |
| `:dvor` | `dvor_dt` | s-2 | `tendency_of_atmosphere_relative_vorticity` | 3D |
| `:ddiv` | `ddiv_dt` | s-2 | `tendency_of_divergence_of_wind` | 3D |
| `:dq` | `dq_dt` | s-1 | `tendency_of_specific_humidity` | 3D |
| **Tracers** |  |  |  |  |
| `:tr1` ... `:tr10` | `tr#` | kg/kg | Passive Tracers 1-10 | 3D |

#### Shallow Water (`:ShallowWater`)

Variables specific to the single-layer shallow water system.

| Symbol | NetCDF Name | Units | Standard Name (`std_name`) | Dimensions |
| :--- | :--- | :--- | :--- | :--- |
| `:h` | `height` | m | `geopotential_height` | 2D |
| `:u` | `ua` | m s-1 | `eastward_wind` | 2D |
| `:v` | `va` | m s-1 | `northward_wind` | 2D |
| `:vor` | `vor` | s-1 | `atmosphere_relative_vorticity` | 2D |
| `:div` | `div` | s-1 | `divergence_of_wind` | 2D |
| `:pv` | `pv` | m-1 s-1 | `atmosphere_potential_vorticity` | 2D |
| `:dh` | `dh_dt` | m s-1 | `tendency_of_geopotential_height` | 2D |

#### Barotropic (`:Barotropic`)

Variables for the purely barotropic vorticity equation.

| Symbol | NetCDF Name | Units | Standard Name (`std_name`) | Dimensions |
| :--- | :--- | :--- | :--- | :--- |
| `:vor` | `vor` | s-1 | `atmosphere_relative_vorticity` | 2D |
| `:u` | `ua` | m s-1 | `eastward_wind` | 2D |
| `:v` | `va` | m s-1 | `northward_wind` | 2D |
| `:ke` | `ke` | m2 s-2 | `specific_kinetic_energy` | 2D |
| `:dvor` | `dvor` | s-2 | `tendency_of_atmosphere_relative_vorticity` | 2D |

---

### Example Configuration Script


```julia
using JGCM

# 1. Physics Configuration
physics_params = Dict{String, Any}(
    
    # Corrections    
    "do_mass_correction"      => true,
    "do_energy_correction"    => true,
    "do_water_correction"     => true,
    "use_virtual_temperature" => true,

    # Grid scale condensation
    "do_Lscale_Cond" => true,
    "L"              => 0.2,

    # PBL fluxes
    "do_Sensible_Heating"    => true,
    "C_H"                    => 0.0044,
    "do_Surface_Evaporation" => true,
    "C_E"                    => 0.0044,
    "do_Implicit_PBL_Scheme" => true,
    "C_D"                    => 0.0044,
    
    "PBL_Top_Mode"  => :PressureLevel,
    "PBL_Top_Value" => 85000.0,

    # Held-Suarez
    "do_HS_Forcing" => true,
    "σ_b"           => 0.7,
    "k_a"           => 1.0/(40.0),
    "k_s"           => 1.0/(4.0),
    "k_f"           => 1.0/(1.0),
    "ΔT_y"          => 60.0, 
    "Δθ_z"          => 10.0
)

# 2. Define Output Paths *Before* Configuration
experiment_name  = "HSt42"
output_path_base = joinpath("exp", experiment_name) 
mkpath(output_path_base)

# 3. Model Configuration
config = Model_Config(

    name = "HS_Moist_T42",
    model_type = :PrimitiveEquation,
    
    # Resolution
    num_fourier = 42, 
    nθ          = 64, 
    nd          = 20, 
    
    # Vertical Coordinate
    vert_coord_option      = "even_sigma", 
    vert_difference_option = "simmons_and_burridge", 
    vert_ref_level_option  = "second_centered_wts",
    
    # Planet Settings
    radius     = 6371.0e3, 
    omega      = 7.292e-5, 
    grav       = 9.80,
    day_to_sec = 86400,
    
    # Time Integration
    Δt         = 600,
    end_time   = 86400 * 4,
    spinup_day = 0.0,
    
    # Numerics
    damping_order = 4, 
    damping_coef  = 1.15741e-4, 
    robert_coef   = 0.04, 
    implicit_coef = 0.5,
    
    # Restart
    is_restart        = false,
    restart_file      = "",
    saving_frequency = 0,    # disable saving restarts

    # Cold start
    initial_condition = :Moist_Spinup,

    # Physics
    moisture_processes = true,
    num_tracers = 1,
    
    # IO
    output_path     = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger          = joinpath(output_path_base, "logger.log"),

    do_plev_output  = false,
    pressure_levels = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    vars_to_output  = [:u, :v, :w, :q, :t, :ps, :shflx, :lhflx, :precip],
    output_interval = 86400,
    
    # Physics
    physics_params = physics_params
)

# 4. Run Simulation
JGCM_Simulate(config)
```
