module Grid_Tracer_Transport_Module

using Base.Threads
using ..Spectral_Spherical_Mesh_Module
using ..Vert_Coordinate_Module
using ..Atmo_Data_Module

export Grid_Tracer_Workspace,
    Advance_Grid_Tracer!,
    Horizontal_Tracer_Advection!,
    Vertical_Tracer_Advection!,
    Restore_Grid_Tracer_Integral!

const DEFAULT_TRACER_CFL = 0.8
const TRACER_ROUNDOFF_FACTOR = 1024.0
const TRACER_THREAD_THRESHOLD = 32_768

"""Preallocated storage for gridpoint finite-volume tracer transport."""
mutable struct Grid_Tracer_Workspace
    q_work::Array{Float64,3}
    q_stage::Array{Float64,3}
    q_half_lambda::Array{Float64,3}
    q_half_phi::Array{Float64,3}
    slope_lambda::Array{Float64,3}
    slope_phi::Array{Float64,3}
    tendency::Array{Float64,3}
    u_face::Array{Float64,3}
    v_face::Array{Float64,3}
    tracer_flux_lambda::Array{Float64,3}
    low_flux_lambda::Array{Float64,3}
    volume_flux_lambda::Array{Float64,3}
    tracer_flux_phi::Array{Float64,3}
    low_flux_phi::Array{Float64,3}
    volume_flux_phi::Array{Float64,3}
    positivity_ratio::Array{Float64,3}
    ppm_slope::Array{Float64,3}
    q_left::Array{Float64,3}
    q_right::Array{Float64,3}
    q6::Array{Float64,3}
    vertical_flux::Array{Float64,3}
    low_vertical_flux::Array{Float64,3}
    zonal_cell_width::Vector{Float64}
    zonal_face_length::Vector{Float64}
    meridional_face_length::Vector{Float64}
    cell_area::Vector{Float64}
    southward_distance::Vector{Float64}
    northward_distance::Vector{Float64}
    horizontal_cfl_by_level::Vector{Float64}
    vertical_cfl_by_level::Vector{Float64}
    chunk_minimum::Vector{Float64}
    chunk_maximum::Vector{Float64}
    pressure_minimum_by_level::Vector{Float64}
    chunk_isfinite::Vector{Bool}
    pressure_isfinite_by_level::Vector{Bool}
    metrics_radius::Float64
    metrics_ready::Bool
end

function Grid_Tracer_Workspace(nλ::Int, nθ::Int, nd::Int)
    cell() = zeros(Float64, nλ, nθ, nd)
    lonface() = zeros(Float64, nλ, nθ, nd)
    latface() = zeros(Float64, nλ, nθ + 1, nd)
    vertface() = zeros(Float64, nλ, nθ, nd + 1)
    return Grid_Tracer_Workspace(
        cell(), cell(), cell(), cell(), cell(), cell(), cell(),
        lonface(), latface(), lonface(), lonface(), lonface(), latface(), latface(), latface(),
        cell(), cell(), cell(), cell(), cell(), vertface(), vertface(),
        zeros(nθ), zeros(nθ), zeros(nθ + 1), zeros(nθ), zeros(nθ), zeros(nθ),
        zeros(nd), zeros(nd), zeros(max(nθ, nd)), zeros(max(nθ, nd)),
        zeros(nd), fill(true, max(nθ, nd)), fill(true, nd), NaN, false,
    )
end

@inline _west(i, nλ) = i == 1 ? nλ : i - 1
@inline _east(i, nλ) = i == nλ ? 1 : i + 1
@inline _pole_shift(i, nλ) = mod1(i + div(nλ, 2), nλ)
@inline _use_tracer_threads(ncells, nchunks) =
    nthreads() > 1 && ncells >= TRACER_THREAD_THRESHOLD && nchunks > 1

"""
Remove only floating-point-scale negative values before they enter another
reconstruction. A larger undershoot remains a hard transport failure.
"""
function _remove_roundoff_undershoots!(tracer, stage)
    qmin = minimum(tracer)
    qmin >= 0.0 && return nothing
    scale = max(maximum(abs, tracer), 1.0)
    tolerance = TRACER_ROUNDOFF_FACTOR * eps(scale)
    qmin >= -tolerance ||
        error("$stage grid tracer transport violated positivity: minimum=$qmin, tolerance=$tolerance")
    @inbounds for index in eachindex(tracer)
        tracer[index] < 0.0 && (tracer[index] = 0.0)
    end
    return nothing
end

function _finish_transport_stage!(workspace, tracer, stage, nchunks)
    tracer_minimum = Inf
    tracer_maximum = 0.0
    values_are_finite = true
    @inbounds for chunk = 1:nchunks
        tracer_minimum = min(tracer_minimum, workspace.chunk_minimum[chunk])
        tracer_maximum = max(tracer_maximum, workspace.chunk_maximum[chunk])
        values_are_finite &= workspace.chunk_isfinite[chunk]
    end
    values_are_finite || error("$stage grid tracer transport produced non-finite values")
    tracer_minimum < 0.0 && _remove_roundoff_undershoots!(tracer, stage)
    return tracer_maximum
end

@inline function _limited_slope(center, backward, forward)
    raw = 0.5 * (forward - backward)
    lower = min(backward, center, forward)
    upper = max(backward, center, forward)
    return copysign(min(abs(raw), 2.0 * (center - lower), 2.0 * (upper - center)), raw)
end

function _prepare_tracer_metrics!(workspace, mesh)
    nθ = length(mesh.θc)
    size(workspace.u_face) == (mesh.nλ, mesh.nθ, mesh.nd) ||
        throw(DimensionMismatch("grid tracer workspace does not match the mesh"))
    workspace.metrics_ready && workspace.metrics_radius == mesh.radius && return nothing
    radius = mesh.radius
    Δλ = 2π / mesh.nλ
    @inbounds for j = 1:nθ
        workspace.zonal_cell_width[j] = radius * mesh.cosθ[j] * Δλ
        workspace.zonal_face_length[j] = radius * (mesh.θe[j+1] - mesh.θe[j])
        workspace.cell_area[j] =
            radius^2 * Δλ * (sin(mesh.θe[j+1]) - sin(mesh.θe[j]))
        workspace.southward_distance[j] = j == 1 ?
            radius * 2(mesh.θc[1] + π / 2) : radius * (mesh.θc[j] - mesh.θc[j-1])
        workspace.northward_distance[j] = j == nθ ?
            radius * 2(π / 2 - mesh.θc[nθ]) : radius * (mesh.θc[j+1] - mesh.θc[j])
    end
    @inbounds for h = 1:nθ+1
        workspace.meridional_face_length[h] = radius * cos(mesh.θe[h]) * Δλ
    end
    workspace.metrics_radius = radius
    workspace.metrics_ready = true
    return nothing
end

function _prepare_horizontal_faces_level!(workspace, u, v, k)
    nλ, nθ, _ = size(u)
    uf, vf = workspace.u_face, workspace.v_face
    fvλ, fvφ = workspace.volume_flux_lambda, workspace.volume_flux_phi
    @inbounds for j = 1:nθ, i = 1:nλ
        uf[i, j, k] = 0.5 * (u[_west(i, nλ), j, k] + u[i, j, k])
        fvλ[i, j, k] = uf[i, j, k] * workspace.zonal_face_length[j]
    end
    @inbounds for i = 1:nλ
        vf[i, 1, k] = 0.0
        vf[i, nθ+1, k] = 0.0
        fvφ[i, 1, k] = 0.0
        fvφ[i, nθ+1, k] = 0.0
    end
    @inbounds for j = 2:nθ, i = 1:nλ
        vf[i, j, k] = 0.5 * (v[i, j-1, k] + v[i, j, k])
        fvφ[i, j, k] = vf[i, j, k] * workspace.meridional_face_length[j]
    end
    return nothing
end

function _prepare_transport_level!(workspace, tracer, u, v, Δp, M, dt, k)
    nλ, nθ, _ = size(u)
    _prepare_horizontal_faces_level!(workspace, u, v, k)
    fvλ, fvφ = workspace.volume_flux_lambda, workspace.volume_flux_phi

    horizontal_cfl = 0.0
    vertical_cfl = 0.0
    tracer_minimum = Inf
    tracer_maximum = 0.0
    pressure_minimum = Inf
    tracer_isfinite = true
    pressure_isfinite = true
    @inbounds for j = 1:nθ, i = 1:nλ
        ie = _east(i, nλ)
        west_volume = fvλ[i, j, k]
        east_volume = fvλ[ie, j, k]
        south_volume = fvφ[i, j, k]
        north_volume = fvφ[i, j+1, k]
        incoming_volume = max(west_volume, 0.0) + max(-east_volume, 0.0) +
                          max(south_volume, 0.0) + max(-north_volume, 0.0)
        horizontal_cfl =
            max(horizontal_cfl, dt * incoming_volume / workspace.cell_area[j])

        distance = v[i, j, k] >= 0.0 ?
            workspace.southward_distance[j] : workspace.northward_distance[j]
        horizontal_cfl = max(horizontal_cfl, dt * abs(v[i, j, k]) / distance)

        layer_pressure = Δp[i, j, k]
        incoming_mass = max(M[i, j, k], 0.0) + max(-M[i, j, k+1], 0.0)
        vertical_cfl = max(vertical_cfl, dt * incoming_mass / layer_pressure)

        tracer_value = tracer[i, j, k]
        if isfinite(tracer_value)
            tracer_minimum = min(tracer_minimum, tracer_value)
            tracer_maximum = max(tracer_maximum, abs(tracer_value))
        else
            tracer_isfinite = false
        end
        if isfinite(layer_pressure)
            pressure_minimum = min(pressure_minimum, layer_pressure)
        else
            pressure_isfinite = false
        end
    end
    workspace.horizontal_cfl_by_level[k] = horizontal_cfl
    workspace.vertical_cfl_by_level[k] = vertical_cfl
    workspace.chunk_minimum[k] = tracer_minimum
    workspace.chunk_maximum[k] = tracer_maximum
    workspace.pressure_minimum_by_level[k] = pressure_minimum
    workspace.chunk_isfinite[k] = tracer_isfinite
    workspace.pressure_isfinite_by_level[k] = pressure_isfinite
    return nothing
end

function _prepare_transport!(workspace, mesh, tracer, u, v, Δp, M, dt, cfl_limit)
    _prepare_tracer_metrics!(workspace, mesh)
    _, _, nd = size(tracer)
    if _use_tracer_threads(length(tracer), nd)
        @threads for k = 1:nd
            _prepare_transport_level!(workspace, tracer, u, v, Δp, M, dt, k)
        end
    else
        for k = 1:nd
            _prepare_transport_level!(workspace, tracer, u, v, Δp, M, dt, k)
        end
    end

    horizontal_cfl = 0.0
    vertical_cfl = 0.0
    tracer_minimum = Inf
    tracer_maximum = 0.0
    pressure_minimum = Inf
    tracer_isfinite = true
    pressure_isfinite = true
    @inbounds for k = 1:nd
        horizontal_cfl = max(horizontal_cfl, workspace.horizontal_cfl_by_level[k])
        vertical_cfl = max(vertical_cfl, workspace.vertical_cfl_by_level[k])
        tracer_minimum = min(tracer_minimum, workspace.chunk_minimum[k])
        tracer_maximum = max(tracer_maximum, workspace.chunk_maximum[k])
        pressure_minimum = min(pressure_minimum, workspace.pressure_minimum_by_level[k])
        tracer_isfinite &= workspace.chunk_isfinite[k]
        pressure_isfinite &= workspace.pressure_isfinite_by_level[k]
    end
    tracer_isfinite && tracer_minimum >= 0.0 ||
        throw(DomainError(tracer_minimum, "input tracer must be finite and non-negative"))
    pressure_isfinite && pressure_minimum > 0.0 || throw(
        DomainError(
            pressure_minimum,
            "layer pressure thickness must be finite and positive",
        ),
    )
    nh = max(1, ceil(Int, horizontal_cfl / cfl_limit))
    nv = max(1, ceil(Int, vertical_cfl / cfl_limit))
    return nh, nv, tracer_maximum
end

"""
Advance a scalar with spherical multidimensional Van Leer transport followed by
nonuniform monotone PPM vertical transport. Zonal Van Leer and vertical PPM
can integrate across every fully swept cell, as in Isca/FMS. The current
advective-form update nevertheless subcycles to enforce the incoming-CFL bound
required by its positivity-preserving low-order baseline. Winds, mass fluxes,
and layer thicknesses are held fixed during the transport step.
"""
function Advance_Grid_Tracer!(
    workspace::Grid_Tracer_Workspace,
    mesh::Spectral_Spherical_Mesh,
    tracer_out::Array{Float64,3},
    tracer_in::Array{Float64,3},
    u::Array{Float64,3},
    v::Array{Float64,3},
    Δp::Array{Float64,3},
    M::Array{Float64,3},
    dt::Real;
    cfl_limit::Float64 = DEFAULT_TRACER_CFL,
)
    0.0 < cfl_limit <= 1.0 || throw(ArgumentError("tracer CFL limit must lie in (0, 1]"))
    dt > 0.0 || throw(ArgumentError("tracer timestep must be positive"))
    size(tracer_in) == size(tracer_out) == size(u) == size(v) == size(Δp) ||
        throw(DimensionMismatch("cell-centered tracer transport fields must have equal sizes"))
    size(M) == (size(tracer_in, 1), size(tracer_in, 2), size(tracer_in, 3) + 1) ||
        throw(DimensionMismatch("vertical mass flux must have nd+1 interfaces"))
    nh, nv, tracer_scale =
        _prepare_transport!(workspace, mesh, tracer_in, u, v, Δp, M, dt, cfl_limit)
    copyto!(workspace.q_work, tracer_in)

    dt_h = dt / nh
    for _ = 1:nh
        tolerance = TRACER_ROUNDOFF_FACTOR * eps(max(tracer_scale, 1.0))
        _horizontal_tracer_advection_prepared!(
            workspace, mesh, workspace.q_stage, workspace.q_work, u, v, dt_h, tolerance,
        )
        tracer_scale =
            _finish_transport_stage!(workspace, workspace.q_stage, "horizontal", size(tracer_in, 3))
        workspace.q_work, workspace.q_stage = workspace.q_stage, workspace.q_work
    end

    dt_v = dt / nv
    for _ = 1:nv
        tolerance = TRACER_ROUNDOFF_FACTOR * eps(max(tracer_scale, 1.0))
        _vertical_tracer_advection!(
            workspace, workspace.q_stage, workspace.q_work, Δp, M, dt_v, tolerance,
        )
        tracer_scale =
            _finish_transport_stage!(workspace, workspace.q_stage, "vertical", size(tracer_in, 2))
        workspace.q_work, workspace.q_stage = workspace.q_stage, workspace.q_work
    end

    copyto!(tracer_out, workspace.q_work)
    return (horizontal_substeps = nh, vertical_substeps = nv)
end

Base.@noinline function _horizontal_positivity_error(
    i, j, k, q, q_low, dt, area,
    west_volume, east_volume, south_volume, north_volume,
)
    meridional_inflow = max(south_volume, 0.0) + max(-north_volume, 0.0)
    zonal_convergence = max(west_volume - east_volume, 0.0)
    zonal_stretching = max(east_volume - west_volume, 0.0)
    incoming = max(west_volume, 0.0) + max(-east_volume, 0.0) + meridional_inflow
    error(
        "horizontal donor-cell update violated positivity at " *
        "(i=$i, j=$j, k=$k): q=$q, q_low=$q_low, dt=$dt, " *
        "west_volume=$west_volume, east_volume=$east_volume, " *
        "south_volume=$south_volume, north_volume=$north_volume, " *
        "incoming_cfl=$(dt * incoming / area), " *
        "zonal_convergence_cfl=$(dt * zonal_convergence / area), " *
        "zonal_stretching_cfl=$(dt * zonal_stretching / area)",
    )
end

function _longitude_slopes!(slope, q)
    nλ, nθ, nd = size(q)
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        slope[i, j, k] = _limited_slope(q[i, j, k], q[_west(i, nλ), j, k], q[_east(i, nλ), j, k])
    end
    return nothing
end

function _latitude_slopes!(slope, q, mesh)
    nλ, nθ, nd = size(q)
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        if j == 1
            qm = q[_pole_shift(i, nλ), 1, k]
            θm = -π - mesh.θc[1]
        else
            qm = q[i, j-1, k]
            θm = mesh.θc[j-1]
        end
        if j == nθ
            qp = q[_pole_shift(i, nλ), nθ, k]
            θp = π - mesh.θc[nθ]
        else
            qp = q[i, j+1, k]
            θp = mesh.θc[j+1]
        end
        qc = q[i, j, k]
        dm = mesh.θc[j] - θm
        dp = θp - mesh.θc[j]
        derivative = (dm * (qp - qc) / dp + dp * (qc - qm) / dm) / (dm + dp)
        raw = derivative * (mesh.θe[j+1] - mesh.θe[j])
        lower = min(qm, qc, qp)
        upper = max(qm, qc, qp)
        slope[i, j, k] = copysign(min(abs(raw), 2(qc - lower), 2(upper - qc)), raw)
    end
    return nothing
end

function _transverse_states!(workspace, mesh, q, u, v, dt)
    qx, qy = workspace.q_half_lambda, workspace.q_half_phi
    nλ, nθ, nd = size(q)
    Δλ = 2π / nλ
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        # Isca's zonal transverse predictor: linearly interpolate q at the
        # semi-Lagrangian departure point. Unlike a one-cell upwind predictor,
        # this remains well-defined when the half-step crosses several cells.
        c = 0.5dt * u[i, j, k] / (mesh.radius * mesh.cosθ[j] * Δλ)
        offset = floor(Int, c)
        fraction = c - offset
        ileft = mod1(i - 1 - offset, nλ)
        iright = _east(ileft, nλ)
        qx[i, j, k] = fraction * q[ileft, j, k] + (1 - fraction) * q[iright, j, k]

        if v[i, j, k] >= 0.0
            qm = j == 1 ? q[_pole_shift(i, nλ), 1, k] : q[i, j-1, k]
            distance = j == 1 ? mesh.radius * 2(mesh.θc[1] + π / 2) :
                       mesh.radius * (mesh.θc[j] - mesh.θc[j-1])
            qy[i, j, k] = q[i, j, k] + 0.5dt * v[i, j, k] * (qm - q[i, j, k]) / distance
        else
            qp = j == nθ ? q[_pole_shift(i, nλ), nθ, k] : q[i, j+1, k]
            distance = j == nθ ? mesh.radius * 2(π / 2 - mesh.θc[nθ]) :
                       mesh.radius * (mesh.θc[j+1] - mesh.θc[j])
            qy[i, j, k] = q[i, j, k] - 0.5dt * v[i, j, k] * (qp - q[i, j, k]) / distance
        end
    end
    return nothing
end

"""
Return low- and high-order averages over the zonal interval swept through a
face. Fully crossed cells contribute their exact cell means; Van Leer is used
only in the fractional terminal donor cell.
"""
@inline function _zonal_swept_averages(q_high, q_low, slope, i, j, k, courant, nλ)
    distance = abs(courant)
    iszero(distance) && return q_low[i, j, k], q_high[i, j, k]
    direction = courant > 0.0 ? -1 : 1
    donor = courant > 0.0 ? _west(i, nλ) : i
    whole = floor(Int, distance)
    fraction = distance - whole
    low_integral = 0.0
    high_integral = 0.0
    @inbounds for _ = 1:whole
        low_integral += q_low[donor, j, k]
        high_integral += q_high[donor, j, k]
        donor = mod1(donor + direction, nλ)
    end
    if fraction > 0.0
        low_value = q_low[donor, j, k]
        high_value = q_high[donor, j, k]
        face_average = courant > 0.0 ?
            high_value + 0.5 * slope[donor, j, k] * (1 - fraction) :
            high_value - 0.5 * slope[donor, j, k] * (1 - fraction)
        low_integral += fraction * low_value
        high_integral += fraction * face_average
    end
    return low_integral / distance, high_integral / distance
end

function _horizontal_level!(workspace, mesh, q_out, q, u, v, dt, tolerance, k)
    nλ, nθ, _ = size(q)
    qx, qy = workspace.q_half_lambda, workspace.q_half_phi
    @inbounds for j = 1:nθ, i = 1:nλ
        c = 0.5dt * u[i, j, k] / workspace.zonal_cell_width[j]
        offset = floor(Int, c)
        fraction = c - offset
        ileft = mod1(i - 1 - offset, nλ)
        iright = _east(ileft, nλ)
        qx[i, j, k] = fraction * q[ileft, j, k] + (1 - fraction) * q[iright, j, k]

        if v[i, j, k] >= 0.0
            qm = j == 1 ? q[_pole_shift(i, nλ), 1, k] : q[i, j-1, k]
            distance = workspace.southward_distance[j]
            qy[i, j, k] =
                q[i, j, k] + 0.5dt * v[i, j, k] * (qm - q[i, j, k]) / distance
        else
            qp = j == nθ ? q[_pole_shift(i, nλ), nθ, k] : q[i, j+1, k]
            distance = workspace.northward_distance[j]
            qy[i, j, k] =
                q[i, j, k] - 0.5dt * v[i, j, k] * (qp - q[i, j, k]) / distance
        end
    end

    @inbounds for j = 1:nθ, i = 1:nλ
        workspace.slope_lambda[i, j, k] = _limited_slope(
            qy[i, j, k], qy[_west(i, nλ), j, k], qy[_east(i, nλ), j, k],
        )

        if j == 1
            qm = qx[_pole_shift(i, nλ), 1, k]
            θm = -π - mesh.θc[1]
        else
            qm = qx[i, j-1, k]
            θm = mesh.θc[j-1]
        end
        if j == nθ
            qp = qx[_pole_shift(i, nλ), nθ, k]
            θp = π - mesh.θc[nθ]
        else
            qp = qx[i, j+1, k]
            θp = mesh.θc[j+1]
        end
        qc = qx[i, j, k]
        dm = mesh.θc[j] - θm
        dp = θp - mesh.θc[j]
        derivative = (dm * (qp - qc) / dp + dp * (qc - qm) / dm) / (dm + dp)
        raw = derivative * (mesh.θe[j+1] - mesh.θe[j])
        lower = min(qm, qc, qp)
        upper = max(qm, qc, qp)
        workspace.slope_phi[i, j, k] =
            copysign(min(abs(raw), 2(qc - lower), 2(upper - qc)), raw)
    end

    uf, vf = workspace.u_face, workspace.v_face
    fqλ, flλ, fvλ = workspace.tracer_flux_lambda,
    workspace.low_flux_lambda,
    workspace.volume_flux_lambda
    fqφ, flφ, fvφ = workspace.tracer_flux_phi,
    workspace.low_flux_phi,
    workspace.volume_flux_phi

    @inbounds for j = 1:nθ, i = 1:nλ
        c = uf[i, j, k] * dt / workspace.zonal_cell_width[j]
        low_average, high_average = _zonal_swept_averages(
            workspace.q_half_phi, q, workspace.slope_lambda, i, j, k, c, nλ,
        )
        fqλ[i, j, k] = fvλ[i, j, k] * high_average
        flλ[i, j, k] = fvλ[i, j, k] * low_average
    end

    @inbounds for i = 1:nλ
        fqφ[i, 1, k] = 0.0
        fqφ[i, nθ+1, k] = 0.0
        flφ[i, 1, k] = 0.0
        flφ[i, nθ+1, k] = 0.0
    end
    @inbounds for h = 2:nθ, i = 1:nλ
        if vf[i, h, k] >= 0.0
            donor = h - 1
            c = vf[i, h, k] * dt / workspace.zonal_face_length[donor]
            qface = workspace.q_half_lambda[i, donor, k] + 0.5 * workspace.slope_phi[i, donor, k] * (1 - c)
        else
            donor = h
            c = -vf[i, h, k] * dt / workspace.zonal_face_length[donor]
            qface = workspace.q_half_lambda[i, donor, k] - 0.5 * workspace.slope_phi[i, donor, k] * (1 - c)
        end
        fqφ[i, h, k] = fvφ[i, h, k] * qface
        flφ[i, h, k] = fvφ[i, h, k] * q[i, donor, k]
    end

    # Flux-corrected transport: the donor-cell update is positive under the
    # combined incoming-CFL bound. Limit only the high-order antidiffusive
    # correction, using one factor per losing cell and one shared face flux.
    ratio = workspace.positivity_ratio
    @inbounds for j = 1:nθ, i = 1:nλ
        ie = _east(i, nλ)
        area = workspace.cell_area[j]
        low_div = flλ[ie, j, k] - flλ[i, j, k] + flφ[i, j+1, k] - flφ[i, j, k]
        volume_div = fvλ[ie, j, k] - fvλ[i, j, k] + fvφ[i, j+1, k] - fvφ[i, j, k]
        q_low = q[i, j, k] + dt * (-low_div + q[i, j, k] * volume_div) / area
        q_low >= -tolerance || _horizontal_positivity_error(
            i, j, k, q[i, j, k], q_low, dt, area,
            fvλ[i, j, k], fvλ[ie, j, k], fvφ[i, j, k], fvφ[i, j+1, k],
        )

        correction_e = fqλ[ie, j, k] - flλ[ie, j, k]
        correction_w = fqλ[i, j, k] - flλ[i, j, k]
        correction_n = fqφ[i, j+1, k] - flφ[i, j+1, k]
        correction_s = fqφ[i, j, k] - flφ[i, j, k]
        loss = max(correction_e, 0.0) + max(-correction_w, 0.0) +
               max(correction_n, 0.0) + max(-correction_s, 0.0)
        ratio[i, j, k] = loss > 0.0 ? min(1.0, max(q_low, 0.0) * area / (dt * loss)) : 1.0
    end

    @inbounds for j = 1:nθ, i = 1:nλ
        correction = fqλ[i, j, k] - flλ[i, j, k]
        losing_i = correction >= 0.0 ? _west(i, nλ) : i
        fqλ[i, j, k] = flλ[i, j, k] + ratio[losing_i, j, k] * correction
    end
    @inbounds for h = 2:nθ, i = 1:nλ
        correction = fqφ[i, h, k] - flφ[i, h, k]
        losing_j = correction >= 0.0 ? h - 1 : h
        fqφ[i, h, k] = flφ[i, h, k] + ratio[i, losing_j, k] * correction
    end

    tracer_minimum = Inf
    tracer_maximum = 0.0
    values_are_finite = true
    @inbounds for j = 1:nθ, i = 1:nλ
        ie = _east(i, nλ)
        area = workspace.cell_area[j]
        tracer_div = fqλ[ie, j, k] - fqλ[i, j, k] + fqφ[i, j+1, k] - fqφ[i, j, k]
        volume_div = fvλ[ie, j, k] - fvλ[i, j, k] + fvφ[i, j+1, k] - fvφ[i, j, k]
        tendency = (-tracer_div + q[i, j, k] * volume_div) / area
        workspace.tendency[i, j, k] = tendency
        tracer_value = q[i, j, k] + dt * tendency
        q_out[i, j, k] = tracer_value
        if isfinite(tracer_value)
            tracer_minimum = min(tracer_minimum, tracer_value)
            tracer_maximum = max(tracer_maximum, abs(tracer_value))
        else
            values_are_finite = false
        end
    end
    workspace.chunk_minimum[k] = tracer_minimum
    workspace.chunk_maximum[k] = tracer_maximum
    workspace.chunk_isfinite[k] = values_are_finite
    return nothing
end

function _prepare_horizontal_faces!(workspace, mesh, u, v)
    _prepare_tracer_metrics!(workspace, mesh)
    _, _, nd = size(u)
    if _use_tracer_threads(length(u), nd)
        @threads for k = 1:nd
            _prepare_horizontal_faces_level!(workspace, u, v, k)
        end
    else
        for k = 1:nd
            _prepare_horizontal_faces_level!(workspace, u, v, k)
        end
    end
    return nothing
end

function _horizontal_tracer_advection_prepared!(workspace, mesh, q_out, q, u, v, dt, tolerance)
    _, _, nd = size(q)
    if _use_tracer_threads(length(q), nd)
        @threads for k = 1:nd
            _horizontal_level!(workspace, mesh, q_out, q, u, v, dt, tolerance, k)
        end
    else
        for k = 1:nd
            _horizontal_level!(workspace, mesh, q_out, q, u, v, dt, tolerance, k)
        end
    end
    return nothing
end

"""One multidimensional spherical Van Leer step with swept-cell zonal fluxes."""
function Horizontal_Tracer_Advection!(workspace, mesh, q_out, q, u, v, dt)
    _prepare_horizontal_faces!(workspace, mesh, u, v)
    tolerance = TRACER_ROUNDOFF_FACTOR * eps(max(maximum(abs, q), 1.0))
    _horizontal_tracer_advection_prepared!(workspace, mesh, q_out, q, u, v, dt, tolerance)
    return nothing
end

Base.@noinline function _vertical_positivity_error(
    i, j, k, q, q_low, dt, Δp, mass_top, mass_bottom,
)
    incoming = max(mass_top, 0.0) + max(-mass_bottom, 0.0)
    error(
        "vertical donor-cell update violated positivity at " *
        "(i=$i, j=$j, k=$k): q=$q, q_low=$q_low, dt=$dt, Δp=$Δp, " *
        "mass_top=$mass_top, mass_bottom=$mass_bottom, " *
        "incoming_cfl=$(dt * incoming / Δp)",
    )
end

function _ppm_reconstruct_column!(workspace, q, dz, i, j)
    nd = size(q, 3)
    s, ql, qr, q6 = workspace.ppm_slope, workspace.q_left, workspace.q_right, workspace.q6
    @inbounds begin
        s[i, j, 1] = 0.0
        s[i, j, nd] = 0.0
        for k = 2:nd-1
            gm = (q[i, j, k] - q[i, j, k-1]) / (dz[i, j, k] + dz[i, j, k-1])
            gp = (q[i, j, k+1] - q[i, j, k]) / (dz[i, j, k+1] + dz[i, j, k])
            raw = (gp * (2dz[i, j, k-1] + dz[i, j, k]) +
                   gm * (2dz[i, j, k+1] + dz[i, j, k])) * dz[i, j, k] /
                  (dz[i, j, k-1] + dz[i, j, k] + dz[i, j, k+1])
            lower = min(q[i, j, k-1], q[i, j, k], q[i, j, k+1])
            upper = max(q[i, j, k-1], q[i, j, k], q[i, j, k+1])
            s[i, j, k] = copysign(
                min(abs(raw), 2(q[i, j, k] - lower), 2(upper - q[i, j, k])), raw,
            )
        end
        for k = 1:nd
            ql[i, j, k] = q[i, j, k] - 0.5s[i, j, k]
            qr[i, j, k] = q[i, j, k] + 0.5s[i, j, k]
        end
        for h = 3:nd-1
            km, kp = h - 1, h
            hmm, hm = dz[i, j, km-1], dz[i, j, km]
            hp, hpp = dz[i, j, kp], dz[i, j, kp+1]
            inv_pair = 1 / (hm + hp)
            inv_four = 1 / (hmm + hm + hp + hpp)
            x = (hmm + hm) / (2hm + hp) - (hp + hpp) / (hm + 2hp)
            y = 2hm * hp
            w1 = hm * inv_pair + x * y * inv_pair * inv_four
            w2 = hm * (hmm + hm) / (2hm + hp) * inv_four
            w3 = hp * (hp + hpp) / (hm + 2hp) * inv_four
            edge = q[i, j, km] + w1 * (q[i, j, kp] - q[i, j, km]) -
                   w2 * s[i, j, kp] + w3 * s[i, j, km]
            qr[i, j, km] = edge
            ql[i, j, kp] = edge
        end
        for k = 1:nd
            qc = q[i, j, k]
            left, right = ql[i, j, k], qr[i, j, k]
            if (right - qc) * (qc - left) <= 0.0
                left = qc
                right = qc
            elseif k != 1 && k != nd
                span = right - left
                curvature = span * (qc - 0.5 * (right + left))
                bound = span^2 / 6
                curvature > bound && (left = 3qc - 2right)
                curvature < -bound && (right = 3qc - 2left)
            end
            ql[i, j, k] = left
            qr[i, j, k] = right
            q6[i, j, k] = 6(qc - 0.5 * (left + right))
        end
    end
    return nothing
end

function _ppm_reconstruction!(workspace, q, dz)
    nλ, nθ, _ = size(q)
    @inbounds for j = 1:nθ, i = 1:nλ
        _ppm_reconstruct_column!(workspace, q, dz, i, j)
    end
    return nothing
end

"""
Return low- and high-order averages over the pressure mass swept through a
vertical interface. Complete donor layers use their exact mass-weighted means;
the terminal partial layer uses its monotone PPM parabola.
"""
@inline function _vertical_swept_averages(workspace, q, Δp, i, j, h, swept_mass)
    total_mass = abs(swept_mass)
    if iszero(total_mass)
        value = q[i, j, clamp(h, 1, size(q, 3))]
        return value, value
    end
    nd = size(q, 3)
    positive = swept_mass > 0.0
    donor = positive ? h - 1 : h
    direction = positive ? -1 : 1
    remaining = total_mass
    low_integral = 0.0
    high_integral = 0.0
    tolerance = 32eps(max(total_mass, 1.0))

    @inbounds while 1 <= donor <= nd && remaining >= Δp[i, j, donor]
        layer_mass = Δp[i, j, donor]
        value = q[i, j, donor]
        low_integral += layer_mass * value
        high_integral += layer_mass * value
        remaining -= layer_mass
        remaining <= tolerance && return low_integral / total_mass, high_integral / total_mass
        donor += direction
    end
    1 <= donor <= nd || error(
        "vertical swept mass crosses the model boundary at interface $h: swept_mass=$swept_mass",
    )

    fraction = remaining / Δp[i, j, donor]
    fraction <= 1.0 + tolerance || error("invalid terminal vertical swept fraction: $fraction")
    fraction = min(fraction, 1.0)
    value = q[i, j, donor]
    if positive
        edge_average = workspace.q_right[i, j, donor] - 0.5fraction * (
            workspace.q_right[i, j, donor] - workspace.q_left[i, j, donor] -
            (1 - 2fraction / 3) * workspace.q6[i, j, donor]
        )
    else
        edge_average = workspace.q_left[i, j, donor] + 0.5fraction * (
            workspace.q_right[i, j, donor] - workspace.q_left[i, j, donor] +
            (1 - 2fraction / 3) * workspace.q6[i, j, donor]
        )
    end
    low_integral += remaining * value
    high_integral += remaining * edge_average
    return low_integral / total_mass, high_integral / total_mass
end

function _vertical_column!(workspace, q_out, q, Δp, M, dt, tolerance, i, j)
    nd = size(q, 3)
    _ppm_reconstruct_column!(workspace, q, Δp, i, j)

    flux, low_flux = workspace.vertical_flux, workspace.low_vertical_flux
    @inbounds begin
        flux[i, j, 1] = 0.0
        flux[i, j, nd+1] = 0.0
        low_flux[i, j, 1] = 0.0
        low_flux[i, j, nd+1] = 0.0
    end
    @inbounds for h = 2:nd
        low_average, high_average = _vertical_swept_averages(
            workspace, q, Δp, i, j, h, dt * M[i, j, h],
        )
        flux[i, j, h] = M[i, j, h] * high_average
        low_flux[i, j, h] = M[i, j, h] * low_average
    end

    ratio = workspace.positivity_ratio
    @inbounds for k = 1:nd
        low_tendency = -(low_flux[i, j, k+1] - low_flux[i, j, k] -
                         q[i, j, k] * (M[i, j, k+1] - M[i, j, k])) / Δp[i, j, k]
        q_low = q[i, j, k] + dt * low_tendency
        q_low >= -tolerance || _vertical_positivity_error(
            i, j, k, q[i, j, k], q_low, dt, Δp[i, j, k],
            M[i, j, k], M[i, j, k+1],
        )
        correction_bottom = flux[i, j, k+1] - low_flux[i, j, k+1]
        correction_top = flux[i, j, k] - low_flux[i, j, k]
        loss = max(correction_bottom, 0.0) + max(-correction_top, 0.0)
        ratio[i, j, k] = loss > 0.0 ?
            min(1.0, max(q_low, 0.0) * Δp[i, j, k] / (dt * loss)) : 1.0
    end
    @inbounds for h = 2:nd
        correction = flux[i, j, h] - low_flux[i, j, h]
        losing_k = correction >= 0.0 ? h - 1 : h
        flux[i, j, h] = low_flux[i, j, h] + ratio[i, j, losing_k] * correction
    end
    tracer_minimum = Inf
    tracer_maximum = 0.0
    values_are_finite = true
    @inbounds for k = 1:nd
        tendency = -(flux[i, j, k+1] - flux[i, j, k] -
                     q[i, j, k] * (M[i, j, k+1] - M[i, j, k])) / Δp[i, j, k]
        workspace.tendency[i, j, k] = tendency
        tracer_value = q[i, j, k] + dt * tendency
        q_out[i, j, k] = tracer_value
        if isfinite(tracer_value)
            tracer_minimum = min(tracer_minimum, tracer_value)
            tracer_maximum = max(tracer_maximum, abs(tracer_value))
        else
            values_are_finite = false
        end
    end
    return tracer_minimum, tracer_maximum, values_are_finite
end

function _vertical_latitude!(workspace, q_out, q, Δp, M, dt, tolerance, j)
    nλ = size(q, 1)
    tracer_minimum = Inf
    tracer_maximum = 0.0
    values_are_finite = true
    @inbounds for i = 1:nλ
        column_minimum, column_maximum, column_isfinite =
            _vertical_column!(workspace, q_out, q, Δp, M, dt, tolerance, i, j)
        tracer_minimum = min(tracer_minimum, column_minimum)
        tracer_maximum = max(tracer_maximum, column_maximum)
        values_are_finite &= column_isfinite
    end
    workspace.chunk_minimum[j] = tracer_minimum
    workspace.chunk_maximum[j] = tracer_maximum
    workspace.chunk_isfinite[j] = values_are_finite
    return nothing
end

function _vertical_tracer_advection!(workspace, q_out, q, Δp, M, dt, tolerance)
    _, nθ, _ = size(q)
    if _use_tracer_threads(length(q), nθ)
        @threads for j = 1:nθ
            _vertical_latitude!(workspace, q_out, q, Δp, M, dt, tolerance, j)
        end
    else
        for j = 1:nθ
            _vertical_latitude!(workspace, q_out, q, Δp, M, dt, tolerance, j)
        end
    end
    return nothing
end

"""One monotone nonuniform-PPM step with multi-layer swept-mass fluxes."""
function Vertical_Tracer_Advection!(workspace, q_out, q, Δp, M, dt)
    tolerance = TRACER_ROUNDOFF_FACTOR * eps(max(maximum(abs, q), 1.0))
    _vertical_tracer_advection!(workspace, q_out, q, Δp, M, dt, tolerance)
    return nothing
end

"""Restore a non-negative grid tracer to an explicitly supplied mass integral."""
function Restore_Grid_Tracer_Integral!(
    tracer,
    grid_ps,
    target_integral,
    vert_coord::Vert_Coordinate,
    mesh::Spectral_Spherical_Mesh,
    atmo_data::Atmo_Data,
)
    isfinite(target_integral) && target_integral >= 0.0 ||
        throw(DomainError(target_integral, "target tracer integral must be finite and non-negative"))
    all(isfinite, tracer) && minimum(tracer) >= -128eps(max(maximum(abs, tracer), 1.0)) ||
        throw(DomainError(minimum(tracer), "tracer must be finite and non-negative before restoration"))
    @. tracer = max(tracer, 0.0)
    if iszero(target_integral)
        fill!(tracer, 0.0)
        return nothing
    end
    current = Mass_Weighted_Global_Integral(vert_coord, mesh, atmo_data, tracer, grid_ps)
    isfinite(current) && current > 0.0 ||
        error("cannot restore a positive tracer integral from a non-positive state")
    tracer .*= target_integral / current
    return nothing
end

end
