using JGCM

const LRF_ARTIFACT = get(
    ENV, "JGCM_LRF_FILE",
    "/home/garywu/undergrad_proposal/LRF/data/latitude_lrf.jld2",
)
const ENABLE_LRF = get(ENV, "JGCM_SMOKE_LRF", "1") == "1"
const RESTART_FILE = get(ENV, "JGCM_SMOKE_RESTART", "")
const OUTPUT_DIRECTORY = get(
    ENV,
    "JGCM_SMOKE_OUTPUT",
    "/data92/garywu/undergrad_proposal/lrf_latitude_smoke_5day",
)

if ENABLE_LRF
    isfile(LRF_ARTIFACT) || error("LRF artifact does not exist: $LRF_ARTIFACT")
end

physics_params = Dict{String,Any}(
    "do_mass_correction" => true,
    "do_energy_correction" => true,
    "do_water_correction" => true,
    "use_virtual_temperature" => true,
    "do_Betts_Miller" => false,
    "initial_humidity_floor" => 0.0,
    "do_Lscale_Cond" => true,
    "condensation_heating_fraction" => 1.0,
    "do_LRF" => ENABLE_LRF,
    "LRF_file" => LRF_ARTIFACT,
    "do_Sensible_Heating" => true,
    "C_H" => 0.0044,
    "do_Surface_Evaporation" => true,
    "C_E" => 0.0044,
    "do_Implicit_PBL_Scheme" => true,
    "C_D" => 0.0044,
    "PBL_Top_Mode" => :PressureLevel,
    "PBL_Top_Value" => 85000.0,
    "do_HS_Forcing" => true,
    "σ_b" => 0.7,
    "k_a" => 1.0 / 40.0,
    "k_s" => 1.0 / 4.0,
    "k_f" => 1.0,
    "T_equator" => 294.0,
    "T_stratosphere" => 200.0,
    "ΔT_y" => 65.0,
    "Δθ_z" => 10.0,
)

config = Model_Config(
    name = "latitude_LRF_smoke",
    institution = "Group of Chaos and Predictability, Department of Atmospheric Sciences, National Taiwan University",
    model_type = :PrimitiveEquation,
    num_fourier = 42,
    nθ = 64,
    nd = 20,
    vert_coord_option = "even_sigma",
    vert_difference_option = "simmons_and_burridge",
    vert_ref_level_option = "second_centered_wts",
    radius = 6371.0e3,
    omega = 7.292e-5,
    grav = 9.80,
    day_to_sec = 86400,
    Δt = 600,
    end_time = 5 * 86400,
    spinup_day = 0.0,
    damping_order = 4,
    damping_coef = 1.15741e-4,
    robert_coef = 0.04,
    implicit_coef = 0.5,
    is_restart = !isempty(RESTART_FILE),
    restart_file = RESTART_FILE,
    saving_frequency = 0,
    initial_condition = :Moist_Spinup,
    moisture_processes = true,
    output_path = OUTPUT_DIRECTORY,
    output_filename = joinpath(OUTPUT_DIRECTORY, "output.nc"),
    logger = joinpath(OUTPUT_DIRECTORY, "logger.log"),
    pressure_levels = Float64[],
    vars_to_output = [:u, :v, :q, :t, :ps, :lrf_dt, :precip],
    output_interval = 86400,
    do_plev_output = false,
    physics_params = physics_params,
)

JGCM_Simulate(config)
