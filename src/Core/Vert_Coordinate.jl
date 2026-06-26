module Vert_Coordinate_Module

import Base
using ..Spectral_Spherical_Mesh_Module
using ..Atmo_Data_Module

export Vert_Coordinate, Compute_Vert_Coord, Vert_Advection!, Mass_Weighted_Global_Integral

"""
There are nd levels (nd+1 interfaces)
a and b should be dimensioned by the number of interfaces = 1 + nd

At these interfaces 
pk = ak + bk * p_surf
here ak is actually ak * p_ref.
where p_ref is a constant reference pressure and p_surf is the instantaneous surface pressure.

"""

mutable struct Vert_Coordinate
    nλ::Int64
    nθ::Int64
    nd::Int64
    vert_coord_option::String
    vert_difference_option::String
    vert_advect_scheme::String

    p_ref::Float64
    zero_top::Bool

    # ak[nd+1] = 0, bk[nd+1] = 1, bk[1] = 0
    ak::Array{Float64,1}
    bk::Array{Float64,1}

    Δak::Array{Float64,1}
    Δbk::Array{Float64,1}

    # memory container
    flux::Array{ComplexF64,3}
    vert_integral::Array{Float64,3}

end



"""
    Vert_Coordinate(
        nλ, nθ, nd,
        vert_coord_option, vert_difference_option, vert_advect_scheme,
        p_ref, zero_top,
        scale_heights, surf_res,
        p_press, p_sigma, exponent
    )

Constructs a `Vert_Coordinate` struct, initializing the vertical coordinate configuration for the atmospheric model.

### Parameters
    - nλ: Longitudinal grid size.
    - nθ: Latitudinal grid size.
    - nd: Number of vertical levels.

    - vert_coord_option: Strategy selector of vertical coordinate (e.g. even_sigma, uneven_sigma, hybrid).
    - vert_difference_option: Selector for vertical finite difference scheme.
    - vert_advect_scheme: Selector for vertical advection algorithm.
    
    - p_ref: Reference pressure (for hybrid coordinate).
    - zero_top: Flag to enforce zero pressure at the top.

    - scale_heights: Parameter for grid stretching/clustering.
    - surf_res: Parameter for surface layer resolution control.

    - p_press: Upper transition threshold for hybrid coordinate.
    - p_sigma: Lower transition threshold for hybrid coordinate.
    - exponent: Power exponent for grid stretching function.

### Returns
    - Vert_Coordinate: A fully initialized data struct.

### Modified
    - nothing

"""
function Vert_Coordinate(
    nλ::Int64,
    nθ::Int64,
    nd::Int64,
    vert_coord_option::String,
    vert_difference_option::String,
    vert_advect_scheme::String,
    p_ref::Float64 = 100000.0,
    zero_top::Bool = true,
    scale_heights::Float64 = 4.0,
    surf_res::Float64 = 1.0,
    p_press::Float64 = 0.1,
    p_sigma::Float64 = 0.3,
    exponent::Float64 = 2.5,
)

    ak, bk = Compute_Vert_Coord(
        nd,
        vert_coord_option,
        p_ref,
        scale_heights,
        surf_res,
        p_press,
        p_sigma,
        exponent,
    )
    Δak, Δbk = ak[2:nd+1] - ak[1:nd], bk[2:nd+1] - bk[1:nd]

    flux = zeros(Float64, nλ, nθ, nd + 1)
    vert_integral = zeros(Float64, nλ, nθ, 1)

    Vert_Coordinate(
        nλ,
        nθ,
        nd,
        vert_coord_option,
        vert_difference_option,
        vert_advect_scheme,
        p_ref,
        zero_top,
        ak,
        bk,
        Δak,
        Δbk,
        flux,
        vert_integral,
    )

end



"""
    Compute_Vert_Coord(
        nd, vert_coord_option,
        p_ref,
        scale_heights, surf_res,
        p_press, p_sigma, exponent
    )

Generates the vertical coordinate coefficients (Ak and Bk) that define the model interfaces.

### Parameters
    - nd: Number of vertical levels.
    - vert_coord_option: Strategy selector of vertical coordinate (e.g. even_sigma, uneven_sigma, hybrid).

    - p_ref: Reference pressure (for hybrid coordinate).

    - scale_heights: Parameter for grid stretching/clustering.
    - surf_res: Parameter for surface layer resolution control.

    - p_press: Upper transition threshold for hybrid coordinate.
    - p_sigma: Lower transition threshold for hybrid coordinate.
    - exponent: Power exponent for grid stretching function.

### Returns
    - a (Array{Float64, 1}): The pressure-dependent interface coefficients (Ak), 
        dimensioned by [nd+1]. Units are Pascals.

    - b (Array{Float64, 1}): The surface-pressure-dependent interface coefficients (Bk), 
        dimensioned by [nd+1]. Units are dimensionless.

### Modified
    - nothing

"""
function Compute_Vert_Coord(
    nd::Int64,
    vert_coord_option::String,
    p_ref::Float64 = 100000.0,
    scale_heights::Float64 = 4.0,
    surf_res::Float64 = 1.0,
    p_press::Float64 = 0.1,
    p_sigma::Float64 = 0.3,
    exponent::Float64 = 2.5,
)

    if (vert_coord_option == "even_sigma")
        a, b = Compute_Even_Sigma(nd)

    elseif (vert_coord_option == "uneven_sigma")
        a, b = Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, true)

    elseif (vert_coord_option == "simmons_and_burridge")
        a, b = Base.invokelatest(Compute_Simmons_Burridge)

    elseif (vert_coord_option == "hybrid")
        a_sigma, b_sigma =
            Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, false)
        b_press, a_press =
            Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, false)
        trans = Transition(nd, b_sigma, p_sigma, p_press)
        a = p_ref * (a_sigma .* trans + a_press .* (1.0 .- trans))
        b = b_sigma .* trans + b_press .* (1.0 .- trans)

    elseif (vert_coord_option == "mcm")
        a, b = Compute_Old_Model_Sigma()

    elseif (vert_coord_option == "v197")
        a, b = Compute_V197_Sigma()

    else
        error("Compute_Vert_Coord ", vert_coord_option, "is not a valid value for option")
    end

    return a, b
end



"""
    Transition(nd, p, p_sigma, p_press)

Computes the vertical transition weights (trans) used to blend sigma and pressure coordinates.

This function generates a weighting profile τ ∈ [0, 1] that dictates the linear combination of
terrain-following (sigma) and fixed (pressure) coordinate definitions. The transition ensures a
smooth morphing from irregular bottom boundaries to flat upper surfaces.

The blending logic is defined by two thresholds:
1. Pure Pressure (τ = 0): When p <= p_press (upper atmosphere). Topography has no effect.
2. Pure Sigma (τ = 1): When p >= p_sigma (lower atmosphere). Grid fully follows terrain.
3. Transition Zone: When p_press < p < p_sigma. Uses a sine-squared interpolation to ensure C1
continuity (smooth derivatives) at the boundaries of the transition zone.

Mathematical formulation:
τ(p) = sin^2( (π/2) * (p - p_press) / (p_sigma - p_press) )

### Parameters
    - nd: Number of vertical levels
    - p: Vertical coordinate values (typically normalized pressure or sigma) at each level.
    - p_sigma: The lower pressure/coordinate threshold. Levels with p >= p_sigma are fully terrain-following (value = 1.0).
    - p_press: The upper pressure/coordinate threshold. Levels with p <= p_press are fully isobaric (value = 0.0).

### Returns
    - trans (Array{Float64, 1})
        The computed weighting factors, dimensioned by [nd+1]. Contains values between 0.0 and 1.0.

### Modified
    - nothing

"""
function Transition(nd::Float64, p::Array{Float64,1}, p_sigma::Float64, p_press::Float64)

    trans = zeros(Float64, nd + 1)

    for k = 1:nd+1
        if (p[k] <= p_press)
            trans[k] = 0.0
        elseif (p[k] >= p_sigma)
            trans[k] = 1.0
        else
            x = p[k] - p_press
            xx = p_sigma - p_press
            trans[k] = (sin(0.5 * pi * x / xx))^2
        end
    end

    return trans
end



"""
    Compute_Even_Sigma(nd)

Generates a uniformly spaced sigma vertical coordinate system.

### Parameters
    - nd: Number of vertical levels.

### Returns
    - a (Array{Float64, 1})
        The pressure-dependent interface coefficients (Ak), filled with zeros. Dimension: [nd+1].

    - b (Array{Float64, 1})
        The surface-pressure-dependent interface coefficients (Bk), linearly spaced from 0.0 to 1.0.
        Dimension: [nd+1].

### Modified
    - nothing

"""
function Compute_Even_Sigma(nd::Int64)

    a = zeros(Float64, nd + 1)
    b = Array(LinRange(0, 1.0, nd + 1))
    return a, b
end



"""
    Compute_Uneven_Sigma(nd, scale_heights, surf_res, exponent, zero_top)

Generates a stretched terrain-following (sigma) vertical coordinate system.

This function creates a non-uniform grid where the vertical levels are clustered to increase
resolution in specific regions (typically the planetary boundary layer). It uses a normalized
coordinate variable (ζ) and a power-law exponent to define a vertical stretching profile.

Mathematical Formulation:
1. Defines a normalized coordinate ζ varying linearly from 1.0 (Top) to 0.0 (Surface).
2. Computes an intermediate height proxy z based on ζ, surface resolution, and an exponent factor.
3. Calculates sigma coefficients (b) using an exponential decay function dependent on 'scale_heights'.

Grid Geometry:
- Top of Atmosphere (k=1): Forces b = 0.0 (if zero_top is true).
- Surface (k=nd+1): Forces b = 1.0.
- The intermediate distribution is controlled by `surf_res` (linear term) and `exponent` (non-linear term).

### Parameters
    - nd: Number of vertical levels.
    - scale_heights: A scaling factor controlling the steepness of the exponential stretching.
    - surf_res: Weighting factor for the linear component of the stretching function. Controls near-surface thickness.
    - exponent: Power parameter applied to the normalized coordinate. Higher values concentrate resolution closer to the surface (or top, depending on the mapping direction).
    - zero_top: If true, explicitly forces the top-most interface coefficient to 0.0, overriding the exponential function result at k=1.

### Returns
    - a (Array{Float64, 1}) The pressure-dependent interface coefficients (Ak), 
        filled with zeros. Dimension: [nd+1].

    - b (Array{Float64, 1})The surface-pressure-dependent interface coefficients (Bk), 
        computed via the stretching function. Dimension: [nd+1].

### Modified
    - nothing

"""
function Compute_Uneven_Sigma(
    nd::Int64,
    scale_heights::Float64,
    surf_res::Float64,
    exponent::Float64,
    zero_top::Bool,
)

    a = zeros(Float64, nd + 1)

    # ζ = 1 - (k-1) / nd
    ζ = Array(LinRange(1.0, 0.0, nd + 1))

    # b = exp( -(surf_res * ζ + (1.0 - surf_res) * (ζ^exponent)) * scale_heights)
    z = -(surf_res * ζ + (1.0 - surf_res) * (ζ .^ exponent))
    b = exp.(z * scale_heights)

    b[nd+1] = 1.0

    if (zero_top)
        b[1] = 0.0
    end
    return a, b

end



"""
    Compute_V197_Sigma()

Loads a pre-defined 18-layer vertical sigma coordinate system (V197 configuration).

Grid Geometry:
- Type: Pure Sigma (Ak = 0 for all levels).
- Resolution: 18 Vertical Layers (nd = 18).
- Distribution: Non-uniform spacing, ranging from 0.0 (top) to 1.0 (surface). The specific values
suggest a resolution distribution tailored for a specific atmospheric regime or model version.

### Parameters
    - nothing

### Returns
    - a (Array{Float64, 1}) The pressure-dependent interface coefficients (Ak), 
        filled with zeros. Dimension: [19].

    - b (Array{Float64, 1}) The surface-pressure-dependent interface coefficients (Bk), 
        containing the hardcoded V197 profile. Dimension: [19].

### Modified
    - nothing

"""
function Compute_V197_Sigma()

    nd = 18

    a = zeros(Float64, nd + 1)
    b = [
        0.0
        0.0089163
        0.0342936
        0.0740741
        0.1262002
        0.1886145
        0.2592592
        0.3360768
        0.4170096
        0.5000000
        0.5829904
        0.6639231
        0.7407407
        0.8113854
        0.8737997
        0.9259259
        0.9657064
        0.9910837
        1.0
    ]


    return a, b
end



"""
    Compute_Old_Model_Sigma()

 Loads a pre-defined 14-layer vertical sigma coordinate system (Old Model configuration).

Grid Geometry:
- Type: Pure Sigma (Ak = 0 for all levels).
- Resolution: 14 Vertical Layers (nd = 14).
- Distribution: Non-uniform spacing, ranging from 0.0 (top) to 1.0 (surface).

### Parameters
    - nothing

### Returns
    - a (Array{Float64, 1}) The pressure-dependent interface coefficients (Ak), 
        filled with zeros. Dimension: [15].

    - b (Array{Float64, 1}) The surface-pressure-dependent interface coefficients (Bk), 
        containing the hardcoded legacy profile. Dimension: [15].

### Modified
    - nothing

"""
function Compute_Old_Model_Sigma()
    nd = 14

    a = zeros(Float64, nd + 1)
    b = [
        0.0
        0.03
        0.0707
        0.1311
        0.2102
        0.3036
        0.4062
        0.5138
        0.6226
        0.7284
        0.8255
        0.9066
        0.9640
        0.9933
        1.0
    ]

    return a, b
end



"""
    Compute_Simmons_Burridge()

Loads the standard Simmons & Burridge (1981) / NCAR CAM 20-layer Hybrid coefficients.

### Parameters
    - nothing

### Returns
    - a (Array{Float64, 1}) The pressure-dependent interface coefficients (Ak),
        containing the hardcoded legacy profile. Dimension: [20].

    - b (Array{Float64, 1}) The surface-pressure-dependent interface coefficients (Bk), 
        containing the hardcoded legacy profile. Dimension: [20].

### Modified
    - nothing

"""
function Compute_Simmons_Burridge()

    ak =
        [
            2.19422,
            4.89520,
            9.88241,
            18.05201,
            29.83724,
            44.62333,
            61.60586,
            78.51243,
            77.31270,
            75.90131,
            74.24086,
            72.28743,
            69.98933,
            67.28574,
            64.10509,
            60.36321,
            55.96111,
            50.78224,
            44.68920,
            37.52196,
            0.0,
        ] .* 100.0 # Convert hPa to Pa

    bk = [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.01505,
        0.03276,
        0.05359,
        0.07810,
        0.10695,
        0.14088,
        0.18080,
        0.22777,
        0.28302,
        0.34801,
        0.42446,
        0.51439,
        1.0,
    ]

    return ak, bk
end



"""
    Vert_Advection(vert_coord, r, dz, w, Δt, vert_advect_scheme, rdt)

Computes the vertical advection tendency for a scalar field using a flux-correction formulation.

Supported Schemes:
- "second_centered": Standard 2nd-order arithmetic mean interpolation (assuming uniform grid).
- "second_centered_wts": 2nd-order interpolation weighted by layer thickness (for variable grids).

### Parameters
    - vert_coord: Vertical grid configuration structure. 
    - r: The scalar field variable to be advected, defined at layer centers (full levels). 
    - dz: Vertical grid spacing or layer thickness (Δz) at each column and level.
    - w: Vertical velocity defined at layer interfaces (half levels).
    - Δt: Time step size.
    - vert_advect_scheme: Selector for the interpolation method used to compute interface values from layer centers.
    - rdt: The tendency array (∂r/∂t) where the result is accumulated.

### Returns
    - nothing

### Modified
    - rdt
    - vert_coord.flux

"""
function Vert_Advection!(
    vert_coord::Vert_Coordinate,
    r::AbstractArray{Float64,3},
    dz::Array{Float64,3},
    w::Array{Float64,3},
    Δt::Int64,
    vert_advect_scheme::String,
    rdt::Array{Float64,3},
)

    nd = vert_coord.nd
    flux = vert_coord.flux

    # no flux boundary condition
    flux[:, :, 1] .= 0.0
    flux[:, :, nd+1] .= 0.0

    # ∂r/∂t = -w * ∂r/∂z
    # -w * ∂r/∂z = -∂(wr)/∂z + r * ∂w/∂z
    # rdt = -w∂r/∂z = -∂wr/∂z + r∂w/∂z = ( [wr]_{k-1/2} -[wr]_{k+1/2})/Δz_k + r_k(w_{k+1/2} - w_{k-1/2})/Δz_k

    # 2nd-order centered scheme assuming variable grid spacing
    # -∂wr/∂z = (wr)_{k-1/2} = w_{k-1/2} * (r_{k-1} + (Δr/Δz)_{k-1/2} * Δz_{k-1})
    if vert_advect_scheme == "second_centered_wts"

        @views begin
            f_in = flux[:, :, 2:nd]
            w_in = w[:, :, 2:nd]
            r_up = r[:, :, 1:nd-1]
            r_dn = r[:, :, 2:nd]
            dz_up = dz[:, :, 1:nd-1]
            dz_dn = dz[:, :, 2:nd]

            @. f_in = w_in * (r_up + (r_dn - r_up) / (dz_up + dz_dn) * dz_up)
        end

        # 2nd-order centered scheme assuming uniform grid spacing
        # -∂wr/∂z = (wr)_{k-1/2} = w_{k-1/2} * (r_{k} + r_{k-1}) * 0.5
    elseif vert_advect_scheme == "second_centered"

        @views begin
            f_in = flux[:, :, 2:nd]
            w_in = w[:, :, 2:nd]
            r_up = r[:, :, 1:nd-1]
            r_dn = r[:, :, 2:nd]

            @. f_in = w_in * (r_up + r_dn) * 0.5
        end

    else
        error("vert_advect_scheme ", vert_advect_scheme, " is not a valid value for option")
    end

    # rdt = ((wr)_{k-1/2} - (wr)_{k+1/2}) / Δz_k + (r_k * (w_{k+1/2} - w_{k-1/2})) / Δz_k
    @views begin
        rdt_k = rdt[:, :, 1:nd]
        flux_k = flux[:, :, 1:nd]
        flux_kp1 = flux[:, :, 2:nd+1]
        r_k = r[:, :, 1:nd]
        w_k = w[:, :, 1:nd]
        w_kp1 = w[:, :, 2:nd+1]
        dz_k = dz[:, :, 1:nd]

        @. rdt_k = (flux_k - flux_kp1 + r_k * (w_kp1 - w_k)) / dz_k
    end

end




"""
    Mass_Weighted_Global_Integral(
        vert_coord, mesh, atmo_data,
        grid_data, grid_ps
    )

Computes the global area-averaged, mass-weighted vertical integral of a scalar field.

### Parameters
    - vert_coord: Vertical grid configuration.
    - mesh: Mesh structure containing geometry information.
    - atmo_data: Physical constants structure.

    - grid_data: The 3D scalar field to be integrated, defined on full model levels.
    - grid_ps: The instantaneous surface pressure field (Pascal).

### Returns
    - mass_weighted_global_integral (Float64)
        The computed global mean value. Units are [units of grid_data] * [kg/m²].

### Modified
    - vert_coord.vert_integral

"""
function Mass_Weighted_Global_Integral(
    vert_coord::Vert_Coordinate,
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
    grid_data::Array{Float64,3},
    grid_ps::Array{Float64,3},
)

    grav = atmo_data.grav
    nd = vert_coord.nd

    Δak, Δbk = vert_coord.Δak, vert_coord.Δbk
    vert_integral = vert_coord.vert_integral

    Δp = similar(grid_ps)

    # Δp_k = Δak + Δbk * p_surface
    # I    = (1/g) * Σ (q_k * Δp_k)
    vert_integral .= 0.0
    for k = 1:nd
        Δp .= Δak[k] .+ Δbk[k] * grid_ps
        vert_integral .+= grid_data[:, :, k] .* Δp[:, :, 1]
    end

    mass_weighted_global_integral = Area_Weighted_Global_Mean(mesh, vert_integral) / grav

    return mass_weighted_global_integral
end

end
