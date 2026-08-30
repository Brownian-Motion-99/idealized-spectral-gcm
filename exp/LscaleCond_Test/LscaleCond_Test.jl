using JGCM

# 1. Resolution (must match Model_Config below)
nθ = 64
nλ = 2 * nθ   # 128

# 2. Build spatially varying latent heat efficiency L[nλ, nθ]
#    L = 0.5 in tropics (|lat| ≤ 30°, i.e. |sinθ| ≤ 0.5), L = 0.0 elsewhere.
sinθ, _ = Compute_Gaussian(nθ)

L_field = zeros(Float64, nλ, nθ)
for j in 1:nθ
    if abs(sinθ[j]) >= 0.5   # sin(30°) = 0.5
        L_field[:, j] .= 0.5
    end
end

# 3. Physics Configuration
physics_params = Dict{String,Any}(

    # Corrections
    "do_mass_correction" => true,
    "do_energy_correction" => true,
    "do_water_correction" => true,
    "use_virtual_temperature" => true,

    # Grid scale condensation — spatially varying efficiency
    "do_Lscale_Cond" => true,
    "condensation_heating_fraction" => L_field,

    # PBL fluxes
    "do_Sensible_Heating" => true,
    "C_H" => 0.0044,
    "do_Surface_Evaporation" => true,
    "C_E" => 0.0044,
    "do_Implicit_PBL_Scheme" => true,
    "C_D" => 0.0044,
    "PBL_Top_Mode" => :PressureLevel,
    "PBL_Top_Value" => 85000.0,

    # Held-Suarez
    "do_HS_Forcing" => true,
    "σ_b" => 0.7,
    "k_a" => 1.0 / (40.0),
    "k_s" => 1.0 / (4.0),
    "k_f" => 1.0 / (1.0),
    "ΔT_y" => 60.0,
    "Δθ_z" => 10.0,
)

# 4. Define Output Paths *Before* Configuration
experiment_name = "LscaleCond_Test"
output_path_base = joinpath("exp", experiment_name)
mkpath(output_path_base)

# 5. Model Configuration
config = Model_Config(
    name = "LscaleCond_Test",
    institution = "Group of Chaos and Predictability, Department of Atmospheric Sciences, National Taiwan University",
    model_type = :PrimitiveEquation,

    # Resolution
    num_fourier = 42,
    nθ = nθ,
    nd = 20,

    # Vertical Coordinate
    vert_coord_option = "even_sigma",
    vert_difference_option = "simmons_and_burridge",
    vert_ref_level_option = "second_centered_wts",

    # Planet Settings
    radius = 6371.0e3,
    omega = 7.292e-5,
    grav = 9.80,
    day_to_sec = 86400,

    # Time Integration — 10-day run, Δt = 600 s
    Δt = 600,
    end_time = 86400 * 10,
    spinup_day = 0.0,

    # Numerics
    damping_order = 4,
    damping_coef = 1.15741e-4,
    robert_coef = 0.04,
    implicit_coef = 0.5,

    # Restart
    is_restart        = false,
    restart_file      = "",
    saving_frequency = 0,

    # Cold start
    initial_condition = :Moist_Spinup,

    # Physics
    moisture_processes = true,
    num_tracers = 1,

    # IO
    output_path = output_path_base,
    output_filename = joinpath(output_path_base, "output.nc"),
    logger          = joinpath(output_path_base, "logger.log"),

    do_plev_output  = true,
    pressure_levels = [100000.0, 92500.0, 85000.0, 70000.0, 50000.0, 30000.0, 20000.0, 10000.0, 5000.0, 1000.0],
    vars_to_output  = [:u, :v, :w, :q, :t, :ps, :shflx, :lhflx, :precip],
    output_interval = 86400,

    # Physics
    physics_params = physics_params,
)

# 6. Run Simulation
JGCM_Simulate(config)
