module PostPhysicsStateSupport

using JGCM

const SD = JGCM.Spectral_Dynamics_Module
const PG = JGCM.Press_And_Geopot_Module
const SI = JGCM.Semi_Implicit_Module
const TI = JGCM.Time_Integrator_Module

function parameters(case::Symbol)
    full = case in (:bm_off, :bm_on)
    Dict{String,Any}(
        "do_mass_correction" => true, "do_energy_correction" => true,
        "do_water_correction" => true, "use_virtual_temperature" => true,
        "initial_humidity_floor" => 0.0,
        "do_Betts_Miller" => case == :bm_on, "bm_tau" => 7200.0,
        "bm_relative_humidity" => 0.7,
        "do_Lscale_Cond" => full, "condensation_heating_fraction" => 1.0,
        "do_Sensible_Heating" => full, "C_H" => 0.0044,
        "do_Surface_Evaporation" => full, "C_E" => 0.0044,
        "do_Implicit_PBL_Scheme" => full || case == :pbl, "C_D" => 0.0044,
        "PBL_Top_Mode" => :PressureLevel, "PBL_Top_Value" => 85000.0,
        "do_HS_Forcing" => full || case in (:hs, :thermal),
        "σ_b" => 0.7, "k_a" => 1 / 40, "k_s" => 1 / 4,
        "k_f" => case == :thermal ? 0.0 : 1.0,
        "T_equator" => 294.0, "T_stratosphere" => 200.0,
        "ΔT_y" => 65.0, "Δθ_z" => 10.0, "do_LRF" => false,
    )
end

function configuration(case; output="/tmp", nf=5, nlat=12, nd=8,
                       dt=600, steps=30)
    Model_Config(
        name="F1_$(case)", model_type=:PrimitiveEquation,
        num_fourier=nf, nθ=nlat, nd=nd, radius=6.371e6,
        omega=7.292e-5, grav=9.8, day_to_sec=86400,
        vert_coord_option="even_sigma", vert_difference_option="simmons_and_burridge",
        vert_ref_level_option="second_centered_wts", Δt=dt, end_time=steps*dt,
        damping_order=4, damping_coef=1.15741e-4, robert_coef=0.04,
        implicit_coef=0.5, initial_condition=:Moist_Spinup,
        moisture_processes=true, saving_frequency=0,
        output_path=output, output_filename=joinpath(output,"output.nc"),
        logger=joinpath(output,"progress.log"), vars_to_output=Symbol[],
        output_interval=43200, physics_params=parameters(case),
    )
end

function with_config(config; kwargs...)
    fields = NamedTuple{fieldnames(Model_Config)}(
        Tuple(getfield(config, key) for key in fieldnames(Model_Config)))
    Model_Config(; merge(fields, (; kwargs...))...)
end

function frame(config; manufactured=true)
    nlon, nlat, nd, nf = 2config.nθ, config.nθ, config.nd, config.num_fourier
    mesh = Spectral_Spherical_Mesh(nf, nf+1, nlon, nlat, nd, config.radius)
    vert = Vert_Coordinate(nlon,nlat,nd,config.vert_coord_option,
        config.vert_difference_option,config.vert_ref_level_option)
    p = config.physics_params
    atmo = Atmo_Data(config.name,nlon,nlat,nd,p["do_mass_correction"],
        p["do_energy_correction"],p["do_water_correction"],
        p["use_virtual_temperature"],mesh.sinθ;
        radius=config.radius,omega=config.omega,grav=config.grav)
    dyn = Dyn_Data(config.name,nf,nf+1,nlon,nlat,nd)
    integrator = TI.Filtered_Leapfrog(config.robert_coef,config.damping_order,
        config.damping_coef,mesh.laplacian_eig,config.implicit_coef,
        config.Δt,true,0,config.end_time)
    semi = SI.Semi_Implicit_Solver(vert,atmo,integrator,1e5,fill(300.0,nd),mesh.wave_numbers)
    bm = JGCM.Driver.betts_miller_state(p,nd,config.Δt,true)
    bm === nothing || (p["BM_state"] = bm)
    Initialize_Atmos_State!(mesh,atmo,dyn,vert,config)
    if manufactured && !config.is_restart
        fill!(dyn.spe_vor_c,0); fill!(dyn.spe_div_c,0)
        for k in 1:nd
            shear = 0.5 + k / nd
            dyn.spe_vor_c[1,2,k] = 2.5e-6 * shear
            dyn.spe_vor_c[2,4,k] = (2e-6 + 0.7e-6im) * shear
            dyn.spe_div_c[3,4,k] = (1e-7 - 0.4e-7im) * shear
        end
        UV_Grid_From_Vor_Div!(mesh,dyn.spe_vor_c,dyn.spe_div_c,dyn.grid_u_c,dyn.grid_v_c)
        Trans_Spherical_To_Grid!(mesh,dyn.spe_vor_c,dyn.grid_vor)
        Trans_Spherical_To_Grid!(mesh,dyn.spe_div_c,dyn.grid_div)
        for key in (:vor,:div,:lnps,:t)
            copyto!(getfield(dyn,Symbol("spe_",key,"_p")),getfield(dyn,Symbol("spe_",key,"_c")))
        end
        for key in (:u,:v,:t,:q,:ps)
            copyto!(getfield(dyn,Symbol("grid_",key,"_p")),getfield(dyn,Symbol("grid_",key,"_c")))
        end
    end
    (;config,mesh,vert,atmo,dyn,integrator,semi)
end

frame(case::Symbol; manufactured=true, kwargs...) =
    frame(configuration(case; kwargs...); manufactured)

function step!(f)
    JGCM.Driver.Step_Dynamics!(f.config,f.mesh,f.atmo,f.dyn,f.integrator,
        f.semi,f.vert,f.config.physics_params)
    f.integrator.init_step && SI.Update_Init_Step!(f.semi)
    f.integrator.time += f.config.Δt
    nothing
end

const PROGNOSTICS = (
    :spe_vor_p,:spe_vor_c,:spe_vor_n,:spe_div_p,:spe_div_c,:spe_div_n,
    :spe_lnps_p,:spe_lnps_c,:spe_lnps_n,:spe_t_p,:spe_t_c,:spe_t_n,
    :grid_u_p,:grid_u_c,:grid_u_n,:grid_v_p,:grid_v_c,:grid_v_n,
    :grid_ps_p,:grid_ps_c,:grid_ps_n,:grid_t_p,:grid_t_c,:grid_t_n,
    :grid_q_p,:grid_q_c,:grid_q_n,
)
snapshot(dyn; fields=PROGNOSTICS) = Dict(key=>copy(getfield(dyn,key)) for key in fields)

function errors!(f, vor, div; level=:c)
    Trans_Spherical_To_Grid!(f.mesh,getfield(f.dyn,Symbol("spe_vor_",level)),vor)
    Trans_Spherical_To_Grid!(f.mesh,getfield(f.dyn,Symbol("spe_div_",level)),div)
    evor = maximum(abs,f.dyn.grid_vor-vor)
    ediv = maximum(abs,f.dyn.grid_div-div)
    (;evor,ediv,svor=maximum(abs,vor),sdiv=maximum(abs,div))
end

end
