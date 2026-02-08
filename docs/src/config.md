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

| Key | Symbol | Default | Physical Meaning |
| :--- | :--- | :--- | :--- |
| `"do_HS_Forcing"` | - | `true` | Master switch for the module. |
| `"σ_b"` |  | `0.7` | Top of the planetary boundary layer (sigma coordinate). |
| `"k_f"` |  | `1.0` day$^{-1}$ | Surface friction damping rate. |
| `"k_a"` |  | `1/40` day$^{-1}$ | Thermal relaxation rate in the free atmosphere. |
| `"k_s"` |  | `1/4` day$^{-1}$ | Thermal relaxation rate at the surface. |
| `"ΔT_y"` |  | `60.0` K | Equator-to-pole temperature difference. |
| `"Δθ_z"` |  | `10.0` K | Vertical potential temperature gradient. |

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
| `end_time` | - | Total simulation duration in seconds. |
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
* **Time:** The model clock resumes exactly where the restart file left off.
* **Filenames:** New restart files will be generated periodically based on `restart_frequency` (in seconds).

---

## Diagnostics & Output

### Standard Output (NetCDF)

The `Output_Manager` writes snapshots to a NetCDF file at `output_interval` (in seconds).

* **Fields:** Defined by `vars_to_output`. See [Available Output Namelist](#available-output-namelist) for available outputs.
* **Vertical Interpolation:** Data is automatically interpolated from model levels to the pressure levels specified in `pressure_levels` (e.g., `[85000.0, 50000.0]`).

### Runtime Logging

The driver prints a status summary to `logger.log` periodically.

* **Monitoring:** It tracks the maximum absolute value of prognostic variables ($U, V, T, P_{surf}$) and their locations.
* **Usage:** Use this log to detect numerical instability (exploding values) early in the run.

---

Based on the `Output_Mappings.jl` file you provided, here is the formatted **Available Output Namelist** section. You can append this to the end of your `configuration.md` or place it within the **Diagnostics & Output** section.

---

### Available Output Namelist

The following tables list the valid symbols available for the `vars_to_output` list in `Model_Config`. The availability depends on the `model_type` selected.

#### Primitive Equation (`:PrimitiveEquation`)

These variables cover full 3D atmospheric states and 2D surface fluxes.

| Symbol | NetCDF Name | Units | Description | Dimensions |
| :--- | :--- | :--- | :--- | :--- |
| **State Variables** |  |  |  |  |
| `:u` | `u` | m/s | Zonal Wind | 3D |
| `:v` | `v` | m/s | Meridional Wind | 3D |
| `:w` | `w` | Pa/s | Vertical Pressure Velocity ($\omega$) | 3D |
| `:t` | `t` | K | Temperature | 3D |
| `:q` | `q` | kg/kg | Specific Humidity | 3D |
| `:vor` | `vor` | 1/s | Relative Vorticity | 3D |
| `:div` | `div` | 1/s | Divergence | 3D |
| `:p` | `p` | Pa | Full Pressure (3D) | 3D |
| `:z` | `z` | m²/s² | Geopotential Height | 3D |
| `:t_eq` | `t_eq` | K | Equilibrium Temperature (HS Forcing) | 3D |
| **Surface & Fluxes** |  |  |  |  |
| `:ps` | `ps` | Pa | Surface Pressure | 2D |
| `:lnps` | `lnps` | numeric | Log Surface Pressure | 2D |
| `:shflx` | `shflx` | W/m² | Sensible Heat Flux | 2D |
| `:lhflx` | `lhflx` | W/m² | Latent Heat Flux | 2D |
| `:precip` | `precip` | mm | Pseudo-adiabatic Precipitation | 2D |
| **Tendencies** |  |  |  |  |
| `:du` | `du_dt` | m/s² | Zonal Wind Tendency | 3D |
| `:dv` | `dv_dt` | m/s² | Meridional Wind Tendency | 3D |
| `:dt` | `dt_dt` | K/s | Temperature Tendency | 3D |
| `:dps` | `dps_dt` | Pa/s | Surface Pressure Tendency | 2D |
| `:dvor` | `dvor_dt` | 1/s² | Relative Vorticity Tendency | 3D |
| `:ddiv` | `ddiv_dt` | 1/s² | Divergence Tendency | 3D |
| `:dq` | `dq_dt` | kg/kg/s | Specific Humidity Tendency | 3D |
| **Tracers** |  |  |  |  |
| `:tr1` ... `:tr10` | `tr#` | kg/kg | Passive Tracers 1-10 | 3D |

#### Shallow Water (`:ShallowWater`)

Variables specific to the single-layer shallow water system.

| Symbol | NetCDF Name | Units | Description |
| :--- | :--- | :--- | :--- |
| `:h` | `h` | m | Geopotential Height |
| `:u` | `u` | m/s | Zonal Wind |
| `:v` | `v` | m/s | Meridional Wind |
| `:vor` | `vor` | 1/s | Relative Vorticity |
| `:div` | `div` | 1/s | Divergence |
| `:pv` | `pv` | s/m | Potential Vorticity |
| `:dh` | `dh_dt` | m/s | Height Tendency |

#### Barotropic (`:Barotropic`)

Variables for the purely barotropic vorticity equation.

| Symbol | NetCDF Name | Units | Description |
| :--- | :--- | :--- | :--- |
| `:vor` | `vor` | 1/s | Relative Vorticity |
| `:u` | `u` | m/s | Zonal Wind |
| `:v` | `v` | m/s | Meridional Wind |
| `:ke` | `ke` | m²/s² | Kinetic Energy |
| `:dvor` | `dvor` | 1/s² | Vorticity Tendency |

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
    restart_frequency = 0,    # disable saving restarts

    # Cold start
    initial_condition = :Moist_Spinup,

    # Physics
    moisture_processes = true,
    num_tracers = 1,
    
    # IO
    output_path     = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger          = joinpath(output_path_base, "logger.log"),

    do_raw_output   = true,
    pressure_levels = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    vars_to_output  = [:u, :v, :w, :q, :t, :ps, :shflx, :lhflx, :precip],
    output_interval = 86400,
    
    # Physics
    physics_params = physics_params
)

# 4. Run Simulation
JGCM_Simulate(config)
```