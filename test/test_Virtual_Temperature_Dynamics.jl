using Test
using JGCM

const Pressure_Geopot_VT = JGCM.Press_And_Geopot_Module
const Spectral_Dynamics_VT = JGCM.Spectral_Dynamics_Module

function four_in_one_state(vert, atmo, grid_tv, p_fields, gradients, u, v, div)
    p_half, Δp, lnp_half, p_full, lnp_full = p_fields
    dλ_ps, dθ_ps = gradients
    nλ, nθ, nd = size(grid_tv)
    M_half = zeros(nλ, nθ, nd + 1)
    w_full = zeros(nλ, nθ, nd)
    δu = zeros(nλ, nθ, nd)
    δv = zeros(nλ, nθ, nd)
    δps = zeros(nλ, nθ, 1)
    δt = zeros(nλ, nθ, nd)
    ps = copy(@view p_half[:, :, end:end])
    Spectral_Dynamics_VT.Four_In_One!(
        vert,
        atmo,
        div,
        u,
        v,
        ps,
        Δp,
        lnp_half,
        lnp_full,
        p_full,
        dλ_ps,
        dθ_ps,
        grid_tv,
        M_half,
        w_full,
        δu,
        δv,
        δps,
        δt,
    )
    return δu, δv, δps, δt, M_half, w_full
end

@testset "Virtual-temperature dynamics consistency" begin
    nλ, nθ, nd = 1, 1, 4
    vert = Vert_Coordinate(
        nλ,
        nθ,
        nd,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )
    moist_atmo = Atmo_Data(
        "moist_operator", nλ, nθ, nd, false, false, false, true, [0.0];
        radius = 6.371e6,
    )
    dry_atmo = Atmo_Data(
        "dry_operator", nλ, nθ, nd, false, false, false, false, [0.0];
        radius = 6.371e6,
    )
    temperature = reshape([225.0, 250.0, 275.0, 300.0], nλ, nθ, nd)
    humidity = reshape([0.001, 0.004, 0.010, 0.020], nλ, nθ, nd)
    virtual_t = similar(temperature)
    dry_t = similar(temperature)
    Pressure_Geopot_VT.Compute_Virtual_Temperature!(
        virtual_t, moist_atmo, temperature, humidity,
    )
    Pressure_Geopot_VT.Compute_Virtual_Temperature!(
        dry_t, dry_atmo, temperature, humidity,
    )
    coefficient = moist_atmo.rvgas / moist_atmo.rdgas - 1.0
    @test virtual_t == temperature .* (1.0 .+ coefficient .* humidity)
    @test dry_t == temperature

    ps = fill(100_000.0, nλ, nθ, 1)
    p_half = zeros(nλ, nθ, nd + 1)
    Δp = zeros(nλ, nθ, nd)
    lnp_half = zeros(nλ, nθ, nd + 1)
    p_full = zeros(nλ, nθ, nd)
    lnp_full = zeros(nλ, nθ, nd)
    Pressure_Geopot_VT.Pressure_Variables!(
        vert, ps, p_half, Δp, lnp_half, p_full, lnp_full,
    )
    p_fields = (p_half, Δp, lnp_half, p_full, lnp_full)
    gradients = (fill(2.5, nλ, nθ, 1), fill(-1.75, nλ, nθ, 1))
    u = fill(7.0, nλ, nθ, nd)
    v = fill(-3.0, nλ, nθ, nd)
    div = fill(1.0e-6, nλ, nθ, nd)
    moist_result = four_in_one_state(
        vert, moist_atmo, virtual_t, p_fields, gradients, u, v, div,
    )
    dry_result = four_in_one_state(
        vert, dry_atmo, dry_t, p_fields, gradients, u, v, div,
    )
    tv_ratio = virtual_t ./ temperature
    @test moist_result[1] ≈ dry_result[1] .* tv_ratio rtol = 2.0e-14
    @test moist_result[2] ≈ dry_result[2] .* tv_ratio rtol = 2.0e-14
    @test moist_result[4] ≈ dry_result[4] .* tv_ratio rtol = 2.0e-14
    @test moist_result[3] == dry_result[3]
    @test moist_result[5] == dry_result[5]
    @test moist_result[6] == dry_result[6]

    zero_humidity = zeros(size(humidity))
    Pressure_Geopot_VT.Compute_Virtual_Temperature!(
        virtual_t, moist_atmo, temperature, zero_humidity,
    )
    @test virtual_t == temperature
    @test_throws DomainError Pressure_Geopot_VT.Compute_Virtual_Temperature!(
        virtual_t, moist_atmo, temperature, fill(1.0, size(humidity)),
    )
end

@testset "Moist hydrostatic pressure-gradient balance" begin
    num_fourier, nθ, nd = 5, 12, 4
    num_spherical = num_fourier + 1
    nλ = 2nθ
    radius = 6.371e6
    mesh = Spectral_Spherical_Mesh(
        num_fourier, num_spherical, nλ, nθ, nd, radius,
    )
    vert = Vert_Coordinate(
        nλ,
        nθ,
        nd,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )
    atmo = Atmo_Data(
        "moist_balance", nλ, nθ, nd, false, false, false, true, mesh.sinθ;
        radius = radius,
    )

    target_tv = 300.0
    humidity_value = 0.02
    coefficient = atmo.rvgas / atmo.rdgas - 1.0
    temperature_value = target_tv / (1.0 + coefficient * humidity_value)
    temperature = fill(temperature_value, nλ, nθ, nd)
    humidity = fill(humidity_value, nλ, nθ, nd)
    virtual_t = similar(temperature)
    Pressure_Geopot_VT.Compute_Virtual_Temperature!(
        virtual_t, atmo, temperature, humidity,
    )

    surface_geopotential = zeros(nλ, nθ, 1)
    for j in 1:nθ
        surface_geopotential[:, j, 1] .= 1000.0 * mesh.sinθ[j]
    end
    grid_lnps = @. log(100_000.0) - surface_geopotential / (atmo.rdgas * target_tv)
    spe_lnps = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, 1)
    Trans_Grid_To_Spherical!(mesh, grid_lnps, spe_lnps)
    Trans_Spherical_To_Grid!(mesh, spe_lnps, grid_lnps)
    ps = exp.(grid_lnps)

    p_half = zeros(nλ, nθ, nd + 1)
    Δp = zeros(nλ, nθ, nd)
    lnp_half = zeros(nλ, nθ, nd + 1)
    p_full = zeros(nλ, nθ, nd)
    lnp_full = zeros(nλ, nθ, nd)
    Pressure_Geopot_VT.Pressure_Variables!(
        vert, ps, p_half, Δp, lnp_half, p_full, lnp_full,
    )
    geopot_full = zeros(nλ, nθ, nd)
    geopot_half = zeros(nλ, nθ, nd + 1)
    Pressure_Geopot_VT.Compute_Geopotential!(
        vert,
        atmo,
        lnp_half,
        lnp_full,
        virtual_t,
        surface_geopotential,
        geopot_full,
        geopot_half,
    )

    dλ_lnps = zeros(nλ, nθ, 1)
    dθ_lnps = zeros(nλ, nθ, 1)
    JGCM.Spectral_Spherical_Mesh_Module.Compute_Gradients!(
        mesh, spe_lnps, dλ_lnps, dθ_lnps,
    )
    dλ_ps = dλ_lnps .* ps
    dθ_ps = dθ_lnps .* ps
    zero_wind = zeros(nλ, nθ, nd)
    pressure_result = four_in_one_state(
        vert,
        atmo,
        virtual_t,
        (p_half, Δp, lnp_half, p_full, lnp_full),
        (dλ_ps, dθ_ps),
        zero_wind,
        zero_wind,
        zero_wind,
    )

    spe_geopot = zeros(ComplexF64, num_fourier + 1, num_spherical + 1, nd)
    dλ_geopot = zeros(nλ, nθ, nd)
    dθ_geopot = zeros(nλ, nθ, nd)
    Trans_Grid_To_Spherical!(mesh, geopot_full, spe_geopot)
    JGCM.Spectral_Spherical_Mesh_Module.Compute_Gradients!(
        mesh, spe_geopot, dλ_geopot, dθ_geopot,
    )
    zonal_residual = pressure_result[1] .- dλ_geopot
    meridional_residual = pressure_result[2] .- dθ_geopot
    gradient_scale = max(maximum(abs, dλ_geopot), maximum(abs, dθ_geopot))
    @test maximum(abs, zonal_residual) <= 1.0e-10 * gradient_scale + 1.0e-15
    @test maximum(abs, meridional_residual) <= 1.0e-10 * gradient_scale + 1.0e-15
end
