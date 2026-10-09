using Test
using JGCM
using JLD2
isdefined(@__MODULE__, :PostPhysicsStateSupport) ||
    include("support/post_physics_state.jl")

@testset "Accepted post-physics vorticity and divergence" begin
    S = PostPhysicsStateSupport
    for case in (:hs,:pbl,:bm_off,:bm_on,:none,:thermal)
        @testset "$case startup and leapfrog" begin
            f = S.frame(case)
            vor, div = similar(f.dyn.grid_vor), similar(f.dyn.grid_div)
            @test maximum(abs,f.dyn.grid_vor) > 1e-7
            @test maximum(abs,f.dyn.grid_div) > 1e-8
            bm_active = false
            for step in 1:30
                S.step!(f)
                e = S.errors!(f,vor,div)
                @test e.evor <= 1e-18 + 5e-13*e.svor
                @test e.ediv <= 1e-18 + 5e-13*e.sdiv
                @test all(isfinite,f.dyn.grid_u_c) && all(isfinite,f.dyn.grid_v_c)
                @test minimum(f.dyn.grid_t_c) > 0 && minimum(f.dyn.grid_ps_c) > 0
                @test minimum(f.dyn.grid_q_c) >= -1e-14
                bm_active |= any(!iszero,f.dyn.grid_bm_t_tendency)
            end
            case == :bm_on && @test bm_active
        end
    end

    @testset "Synchronization uses next spectra and preserves time history" begin
        f = S.frame(:hs)
        S.SD.Spectral_Dynamics!(f.config,f.mesh,f.vert,f.atmo,f.dyn,f.semi;
            advance_time=false)
        S.PG.Pressure_Variables!(f.vert,f.dyn.grid_ps_n,f.dyn.grid_p_half,
            f.dyn.grid_Δp,f.dyn.grid_lnp_half,f.dyn.grid_p_full,f.dyn.grid_lnp_full)
        provisional_u = copy(f.dyn.grid_u_n)
        JGCM.Atmos_Param_Module.Spectral_Physics!(f.config,f.mesh,f.vert,f.atmo,
            f.dyn,f.semi,f.config.physics_params)
        @test maximum(abs,f.dyn.grid_u_n-provisional_u) > 1e-4
        oldfields = filter(key->!endswith(string(key),"_n"),S.PROGNOSTICS)
        before = S.snapshot(f.dyn; fields=oldfields)
        S.SD._synchronize_physics_next!(f.mesh,f.vert,f.atmo,f.dyn)
        vor,div = similar(f.dyn.grid_vor),similar(f.dyn.grid_div)
        e = S.errors!(f,vor,div;level=:n)
        @test e.evor <= 1e-18 + 5e-13*e.svor
        @test e.ediv <= 1e-18 + 5e-13*e.sdiv
        @test all(getfield(f.dyn,key)==value for (key,value) in before)
        Trans_Spherical_To_Grid!(f.mesh,f.dyn.spe_vor_c,vor)
        @test maximum(abs,f.dyn.grid_vor-vor) > 1e-10
    end

    @testset "Rayleigh drag scales resolved curl and divergence" begin
        f = S.frame(:hs)
        d = f.dyn
        copyto!(d.grid_u_n,d.grid_u_c); copyto!(d.grid_v_n,d.grid_v_c)
        fill!(d.grid_ps_n,1e5); fill!(d.grid_t_n,300); fill!(d.grid_q_n,0)
        S.PG.Pressure_Variables!(f.vert,d.grid_ps_n,d.grid_p_half,d.grid_Δp,
            d.grid_lnp_half,d.grid_p_full,d.grid_lnp_full)
        expected_vor, expected_div = copy(d.spe_vor_c), copy(d.spe_div_c)
        for k in 1:f.mesh.nd
            sigma = d.grid_p_full[1,1,k]/1e5
            @test all(d.grid_p_full[:,:,k]./d.grid_p_half[:,:,end] .≈ sigma)
            rate = max(0.0,(sigma-0.7)/0.3)/86400
            factor = 1-600rate
            expected_vor[:,:,k] .*= factor
            expected_div[:,:,k] .*= factor
        end
        JGCM.Atmos_Param_Module.Rayleigh_Friction!(f.atmo,600,86400,
            d.grid_p_half,d.grid_p_full,d.grid_u_n,d.grid_v_n,d.grid_t_n,
            f.config.physics_params)
        before = S.snapshot(d)
        S.SD._synchronize_physics_next!(f.mesh,f.vert,f.atmo,d)
        @test d.spe_vor_n ≈ expected_vor rtol=1e-12 atol=1e-19
        @test d.spe_div_n ≈ expected_div rtol=1e-12 atol=1e-19
        vor,div = similar(d.grid_vor),similar(d.grid_div)
        Trans_Spherical_To_Grid!(f.mesh,expected_vor,vor)
        Trans_Spherical_To_Grid!(f.mesh,expected_div,div)
        @test d.grid_vor ≈ vor rtol=1e-12 atol=1e-18
        @test d.grid_div ≈ div rtol=1e-12 atol=1e-18
        # Reconstruction itself must not mutate any accepted prognostic.
        accepted = S.snapshot(d)
        S.SD._refresh_grid_vor_div!(f.mesh,d.spe_vor_n,d.spe_div_n,d.grid_vor,d.grid_div)
        @test all(getfield(d,key)==value for (key,value) in accepted)
        @test before[:spe_vor_c] == d.spe_vor_c
    end

    @testset "Restart reconstructs diagnostics and preserves trajectory" begin
        mktempdir() do dir
            f = S.frame(:hs;output=dir,steps=5)
            for _ in 1:3; S.step!(f); end
            original = S.snapshot(f.dyn)
            # Model a complete pre-fix checkpoint with stale derived arrays.
            fill!(f.dyn.grid_vor,7e-4); fill!(f.dyn.grid_div,-8e-4)
            manager = Restart_Manager(dir,600)
            Write_Restart_File(manager,f.dyn,f.integrator.time)
            path = joinpath(dir,"restart_t1800.jld2")
            for (key,val) in original; copyto!(getfield(f.dyn,key),val); end
            Trans_Spherical_To_Grid!(f.mesh,f.dyn.spe_vor_c,f.dyn.grid_vor)
            Trans_Spherical_To_Grid!(f.mesh,f.dyn.spe_div_c,f.dyn.grid_div)

            warm = S.with_config(S.configuration(:hs;output=dir,steps=2);
                is_restart=true,restart_file=path)
            loaded = S.frame(warm)
            @test all(getfield(loaded.dyn,key)==val for (key,val) in original)
            vor,div = similar(loaded.dyn.grid_vor),similar(loaded.dyn.grid_div)
            e = S.errors!(loaded,vor,div)
            @test e.evor <= 1e-18 + 5e-13*e.svor
            @test e.ediv <= 1e-18 + 5e-13*e.sdiv

            # Exercise the real driver path, including resumed leapfrog setup.
            run_dir = joinpath(dir,"driver"); mkpath(run_dir)
            driver_config = S.with_config(warm;output_path=run_dir,
                output_filename=joinpath(run_dir,"output.nc"),
                logger=joinpath(run_dir,"progress.log"),saving_frequency=600)
            JGCM_Simulate(driver_config)
            resumed = Dyn_Data("resumed",5,6,24,12,8)
            @test Load_Restart_File!(resumed,
                joinpath(run_dir,"restart","restart_t3000.jld2")) == 3000
            for _ in 1:2; S.step!(f); end
            for key in S.PROGNOSTICS
                @test getfield(resumed,key) ≈ getfield(f.dyn,key) rtol=1e-12 atol=1e-14
            end
            @test resumed.grid_vor ≈ f.dyn.grid_vor rtol=5e-13 atol=1e-18
            @test resumed.grid_div ≈ f.dyn.grid_div rtol=5e-13 atol=1e-18
        end
    end
end
