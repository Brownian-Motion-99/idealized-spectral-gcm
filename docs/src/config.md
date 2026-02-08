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
# 1. Physics
physics_params = Dict{String, Any}(
    "do_mass_correction" => true,
    "do_HS_Forcing"      => true,
    "σ_b"                => 0.7
)

# 2. Config
config = Model_Config(
    name        = "HSt21",
    model_type  = :PrimitiveEquation,
    num_fourier = 21, nθ = 32, nd = 20,
    vert_coord_option = "simmons_and_burridge",
    
    Δt = 600, end_time = 86400 * 10,
    
    is_restart   = false,
    output_path  = "exp/HSt21",
    physics_params = physics_params
)

# 3. Run
JGCM_Simulate(config)

```