using Test
using Random
using JGCM

const Pressure_Geopot = JGCM.Press_And_Geopot_Module
const Vert_Coordinate_Kernels = JGCM.Vert_Coordinate_Module
const Spectral_Mesh_Kernels = JGCM.Spectral_Spherical_Mesh_Module

function reference_pressure_variables(vert, grid_ps)
    nd = vert.nd
    grid_p_half =
        reshape(vert.ak, 1, 1, nd + 1) .+
        reshape(vert.bk, 1, 1, nd + 1) .* grid_ps
    grid_Δp = grid_p_half[:, :, 2:nd+1] - grid_p_half[:, :, 1:nd]
    grid_lnp_half = zeros(size(grid_p_half))
    grid_lnp_full = zeros(size(grid_Δp))
    k_top = vert.zero_top ? 2 : 1

    grid_lnp_half[:, :, k_top:nd+1] .= log.(grid_p_half[:, :, k_top:nd+1])
    grid_lnp_full[:, :, k_top:nd] .=
        grid_lnp_half[:, :, k_top+1:nd+1] .+
        grid_p_half[:, :, k_top:nd] .*
        (grid_lnp_half[:, :, k_top+1:nd+1] - grid_lnp_half[:, :, k_top:nd]) ./
        grid_Δp[:, :, k_top:nd] .- 1.0

    if vert.zero_top
        grid_lnp_half[:, :, 1] .= 0.0
        grid_lnp_full[:, :, 1] .= grid_lnp_half[:, :, 2] .- 1.0
    end

    return grid_p_half, grid_Δp, grid_lnp_half, exp.(grid_lnp_full), grid_lnp_full
end

function reference_geopotential(
    vert,
    atmo,
    grid_lnp_half,
    grid_lnp_full,
    grid_t,
    grid_geopots,
    grid_q,
)
    nd = vert.nd
    grid_geopot_half = zeros(size(grid_lnp_half))
    grid_geopot_full = zeros(size(grid_lnp_full))
    grid_geopot_half[:, :, nd+1] .= grid_geopots[:, :, 1]
    k_top = vert.zero_top ? 2 : 1

    if vert.zero_top
        grid_geopot_half[:, :, 1] .= 0.0
    end

    virtual_t = if atmo.use_virtual_temperature
        grid_t .* (1.0 .+ (atmo.rvgas / atmo.rdgas - 1.0) .* grid_q)
    else
        copy(grid_t)
    end

    for k = nd:-1:k_top
        grid_geopot_half[:, :, k] .=
            grid_geopot_half[:, :, k+1] .+
            atmo.rdgas .* virtual_t[:, :, k] .*
            (grid_lnp_half[:, :, k+1] - grid_lnp_half[:, :, k])
    end
    for k = 1:nd
        grid_geopot_full[:, :, k] .=
            grid_geopot_half[:, :, k+1] .+
            atmo.rdgas .* virtual_t[:, :, k] .*
            (grid_lnp_half[:, :, k+1] - grid_lnp_full[:, :, k])
    end

    return grid_geopot_full, grid_geopot_half, virtual_t
end

function reference_mass_integral(vert, mesh, atmo, grid_data, grid_ps)
    vertical_integral = zeros(size(grid_ps))
    Δp = similar(grid_ps)
    for k = 1:vert.nd
        Δp .= vert.Δak[k] .+ vert.Δbk[k] * grid_ps
        vertical_integral .+= grid_data[:, :, k] .* Δp[:, :, 1]
    end
    area_mean = sum(vertical_integral[:, :, 1] * mesh.wts) / (2mesh.nλ)
    return area_mean / atmo.grav
end

@testset "Allocation-free pressure and geopotential diagnostics" begin
    rng = MersenneTwister(42)
    num_fourier, nθ, nd = 7, 12, 6
    num_spherical = num_fourier + 1
    nλ = 2nθ
    radius = 6.371e6
    mesh = Spectral_Spherical_Mesh(
        num_fourier,
        num_spherical,
        nλ,
        nθ,
        nd,
        radius,
    )
    vert = Vert_Coordinate(
        nλ,
        nθ,
        nd,
        "even_sigma",
        "simmons_and_burridge",
        "second_centered_wts",
    )

    grid_ps = 90_000.0 .+ 15_000.0 .* rand(rng, nλ, nθ, 1)
    grid_t = 220.0 .+ 90.0 .* rand(rng, nλ, nθ, nd)
    grid_q = 0.02 .* rand(rng, nλ, nθ, nd)
    grid_geopots = 5_000.0 .* rand(rng, nλ, nθ, 1)

    p_half_ref, Δp_ref, lnp_half_ref, p_full_ref, lnp_full_ref =
        reference_pressure_variables(vert, grid_ps)
    p_half = similar(p_half_ref)
    Δp = similar(Δp_ref)
    lnp_half = similar(lnp_half_ref)
    p_full = similar(p_full_ref)
    lnp_full = similar(lnp_full_ref)

    Pressure_Geopot.Pressure_Variables!(
        vert,
        grid_ps,
        p_half,
        Δp,
        lnp_half,
        p_full,
        lnp_full,
    )

    @test p_half == p_half_ref
    @test Δp == Δp_ref
    @test lnp_half == lnp_half_ref
    @test lnp_full == lnp_full_ref
    @test p_full == p_full_ref

    pressure_allocations = @allocated Pressure_Geopot.Pressure_Variables!(
        vert,
        grid_ps,
        p_half,
        Δp,
        lnp_half,
        p_full,
        lnp_full,
    )
    @test pressure_allocations == 0

    for use_virtual_temperature in (false, true)
        atmo = Atmo_Data(
            "diagnostic_test",
            nλ,
            nθ,
            nd,
            true,
            true,
            true,
            use_virtual_temperature,
            mesh.sinθ;
            radius = radius,
        )
        geopot_full_ref, geopot_half_ref, virtual_t = reference_geopotential(
            vert,
            atmo,
            lnp_half,
            lnp_full,
            grid_t,
            grid_geopots,
            grid_q,
        )
        geopot_full = similar(geopot_full_ref)
        geopot_half = similar(geopot_half_ref)

        Pressure_Geopot.Compute_Geopotential!(
            vert,
            atmo,
            lnp_half,
            lnp_full,
            grid_t,
            grid_geopots,
            geopot_full,
            geopot_half,
            grid_q,
        )

        @test geopot_full == geopot_full_ref
        @test geopot_half == geopot_half_ref

        # Verify the discrete hydrostatic equations independently at every
        # level, including the virtual-temperature correction when enabled.
        k_top = vert.zero_top ? 2 : 1
        for k = k_top:nd
            @test geopot_half[:, :, k] - geopot_half[:, :, k+1] ≈
                  atmo.rdgas .* virtual_t[:, :, k] .*
                  (lnp_half[:, :, k+1] - lnp_half[:, :, k])
        end
        for k = 1:nd
            @test geopot_full[:, :, k] - geopot_half[:, :, k+1] ≈
                  atmo.rdgas .* virtual_t[:, :, k] .*
                  (lnp_half[:, :, k+1] - lnp_full[:, :, k])
        end

        geopotential_allocations = @allocated Pressure_Geopot.Compute_Geopotential!(
            vert,
            atmo,
            lnp_half,
            lnp_full,
            grid_t,
            grid_geopots,
            geopot_full,
            geopot_half,
            grid_q,
        )
        @test geopotential_allocations == 0
    end
end

@testset "Allocation-free conservation integrals" begin
    rng = MersenneTwister(84)
    num_fourier, nθ, nd = 7, 12, 6
    num_spherical = num_fourier + 1
    nλ = 2nθ
    radius = 6.371e6
    mesh = Spectral_Spherical_Mesh(
        num_fourier,
        num_spherical,
        nλ,
        nθ,
        nd,
        radius,
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
        "integral_test",
        nλ,
        nθ,
        nd,
        true,
        true,
        true,
        true,
        mesh.sinθ;
        radius = radius,
    )

    grid_ps = 90_000.0 .+ 15_000.0 .* rand(rng, nλ, nθ, 1)
    grid_data = randn(rng, nλ, nθ, nd)

    mean_ref = sum(grid_ps[:, :, 1] * mesh.wts) / (2mesh.nλ)
    mean_result = Spectral_Mesh_Kernels.Area_Weighted_Global_Mean(mesh, grid_ps)
    @test mean_result == mean_ref

    integral_ref = reference_mass_integral(vert, mesh, atmo, grid_data, grid_ps)
    integral_result = Vert_Coordinate_Kernels.Mass_Weighted_Global_Integral(
        vert,
        mesh,
        atmo,
        grid_data,
        grid_ps,
    )
    @test integral_result == integral_ref

    # A vertically and horizontally constant field has the analytic integral
    # c * <surface pressure> / g for a sigma coordinate.
    constant_value = 3.25
    constant_field = fill(constant_value, nλ, nθ, nd)
    constant_integral = Vert_Coordinate_Kernels.Mass_Weighted_Global_Integral(
        vert,
        mesh,
        atmo,
        constant_field,
        grid_ps,
    )
    @test constant_integral ≈ constant_value * mean_result / atmo.grav rtol = 2e-15

    # Scalar returns may be boxed by @allocated on Julia 1.8; the former array
    # allocations were much larger, so this bound catches their reintroduction.
    mean_allocations =
        @allocated Spectral_Mesh_Kernels.Area_Weighted_Global_Mean(mesh, grid_ps)
    integral_allocations = @allocated Vert_Coordinate_Kernels.Mass_Weighted_Global_Integral(
        vert,
        mesh,
        atmo,
        grid_data,
        grid_ps,
    )
    @test mean_allocations <= 16
    @test integral_allocations <= 16
end
