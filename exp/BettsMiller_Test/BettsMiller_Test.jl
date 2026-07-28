using JGCM

experiment_name = "BettsMiller_Test"
output_path = joinpath("exp", experiment_name)
mkpath(output_path)

physics_params = Dict{String,Any}(
    "do_mass_correction" => true,
    "do_energy_correction" => true,
    "do_water_correction" => true,
    "use_virtual_temperature" => true,
    "do_Betts_Miller" => true,
    "bm_tau" => 7200.0,
    "bm_relative_humidity" => 0.8,
    # Keeps the spectrally truncated analytic initialization nonnegative.
    "initial_humidity_floor" => 5.0e-5,
    # Phase one deliberately forbids using this direct-state scheme with BM.
    "do_Lscale_Cond" => false,
    "L" => 0.2,
    "do_LRF" => false,
    "do_Sensible_Heating" => true,
    "C_H" => 0.0044,
    "do_Surface_Evaporation" => true,
    "C_E" => 0.0044,
    "do_Implicit_PBL_Scheme" => true,
    "C_D" => 0.0044,
    "PBL_Top_Mode" => :PressureLevel,
    "PBL_Top_Value" => 85000.0,
    "do_HS_Forcing" => true,
    "sigma_b" => 0.7,
    "k_a" => 1.0 / 40.0,
    "k_s" => 1.0 / 4.0,
    "k_f" => 1.0,
    "delta_T_y" => 60.0,
    "delta_theta_z" => 10.0,
)

# HS_Forcing currently uses the legacy Unicode dictionary keys.
physics_params["σ_b"] = physics_params["sigma_b"]
physics_params["ΔT_y"] = physics_params["delta_T_y"]
physics_params["Δθ_z"] = physics_params["delta_theta_z"]

config = Model_Config(
    name = experiment_name,
    institution = "Group of Chaos and Predictability, Department of Atmospheric Sciences, National Taiwan University",
    model_type = :PrimitiveEquation,
    num_fourier = 21,
    nθ = 32,
    nd = 20,
    vert_coord_option = "even_sigma",
    vert_difference_option = "simmons_and_burridge",
    vert_ref_level_option = "second_centered_wts",
    radius = 6371.0e3,
    omega = 7.292e-5,
    grav = 9.80,
    Δt = 600,
    end_time = 86400,
    day_to_sec = 86400,
    damping_order = 4,
    damping_coef = 1.15741e-4,
    robert_coef = 0.04,
    implicit_coef = 0.5,
    initial_condition = :Moist_Spinup,
    moisture_processes = true,
    num_tracers = 1,
    output_path = output_path,
    output_filename = joinpath(output_path, "output.nc"),
    logger = joinpath(output_path, "logger.log"),
    do_plev_output = true,
    vars_to_output = [:u, :v, :q, :t, :ps, :precip, :bm_dt, :bm_dq, :bm_precip],
    output_interval = 3600,
    saving_frequency = 0,
    physics_params = physics_params,
)

JGCM_Simulate(config)
