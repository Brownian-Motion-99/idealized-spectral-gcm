# Getting started

## Install from a clone

Clone the repository, instantiate its Julia environment, and run the tests:

```shell
git clone https://github.com/Brownian-Motion-99/idealized-spectral-gcm.git
cd idealized-spectral-gcm
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

The package entry point is `JGCM_Simulate(config::Model_Config)`. Experiments
are ordinary Julia scripts, so there is no separate namelist parser.

## Run a supplied experiment

The quickest primitive-equation example is the T21 Held--Suarez script:

```shell
julia --project=. exp/HSt21/HS.jl
```

Other supplied examples exercise distinct parts of the model:

| Script | Equation set or feature |
|:---|:---|
| `exp/HSt42/HS.jl` | moist primitive equations, large-scale condensation, surface exchange, and PBL mixing |
| `exp/BettsMiller_Test/BettsMiller_Test.jl` | Betts--Miller convection |
| `exp/LscaleCond_Test/LscaleCond_Test.jl` | spatially varying condensation-heating fraction |
| `exp/Shallow_Water/SW.jl` | shallow-water model |
| `exp/Barotropic/Barotropic.jl` | barotropic vorticity model |

Read an example before running it: output paths and simulation lengths are set
inside each script. In particular, `exp/HSt42/HS.jl` is configured as a long
research run and uses a site-specific absolute output path; adapt both for your
machine.

## Anatomy of an experiment

A primitive-equation experiment has three parts.

First, define parameterized physics:

```julia
using JGCM

physics_params = Dict{String,Any}(
    "do_mass_correction" => true,
    "do_energy_correction" => true,
    "do_water_correction" => true,
    "use_virtual_temperature" => true,

    "do_Betts_Miller" => false,
    "do_Lscale_Cond" => false,
    "do_Sensible_Heating" => false,
    "do_Surface_Evaporation" => false,
    "do_Implicit_PBL_Scheme" => false,

    "do_HS_Forcing" => true,
    "σ_b" => 0.7,
    "k_a" => 1 / 40,
    "k_s" => 1 / 4,
    "k_f" => 1.0,
    "T_equator" => 294.0,
    "T_stratosphere" => 200.0,
    "ΔT_y" => 60.0,
    "Δθ_z" => 10.0,
)
```

Second, construct a complete `Model_Config`:

```julia
run_dir = joinpath("exp", "my_dry_hs")

config = Model_Config(
    name = "my_dry_hs",
    institution = "My institution",
    model_type = :PrimitiveEquation,

    num_fourier = 21,
    nθ = 32,
    nd = 10,
    radius = 6.371e6,
    omega = 7.292e-5,
    grav = 9.8,

    vert_coord_option = "even_sigma",
    vert_difference_option = "simmons_and_burridge",
    vert_ref_level_option = "second_centered_wts",

    Δt = 600,
    end_time = 4 * 86400,
    day_to_sec = 86400,
    damping_order = 4,
    damping_coef = 1.15741e-4,
    robert_coef = 0.04,
    implicit_coef = 0.5,

    initial_condition = :Moist_Spinup,
    moisture_processes = false,
    is_restart = false,
    saving_frequency = 0,

    output_path = run_dir,
    output_filename = joinpath(run_dir, "output.nc"),
    logger = joinpath(run_dir, "logger.log"),
    vars_to_output = [:u, :v, :t, :ps, :t_eq],
    output_interval = 6 * 3600,
    do_plev_output = false,

    physics_params = physics_params,
)
```

Finally, run it:

```julia
JGCM_Simulate(config)
```

The driver creates the configured output directories. NetCDF output is named
from `output_filename`, with a `_t<seconds>` suffix added for the first and any
subsequent chunks. See [Output and restarts](@ref) before enabling a cold start
in a directory that already contains checkpoints.

## Choosing a model

- Use `:PrimitiveEquation` for three-dimensional atmospheric circulation and
  any parameterized physics.
- Use `:Shallow_Water` for a divergent one-layer system with height relaxation
  and momentum damping.
- Use `:Barotropic` for a nondivergent one-layer vorticity system.

The model choice changes valid initial conditions and output symbols. It does
not merely remove vertical levels from the primitive-equation equations.

## Building the documentation locally

Instantiate the documentation environment and run Documenter:

```shell
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The rendered site is written to `docs/build/`.

