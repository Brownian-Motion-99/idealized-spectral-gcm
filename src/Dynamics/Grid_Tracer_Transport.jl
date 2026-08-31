module Grid_Tracer_Transport_Module

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
    )
end

@inline _west(i, nλ) = i == 1 ? nλ : i - 1
@inline _east(i, nλ) = i == nλ ? 1 : i + 1
@inline _pole_shift(i, nλ) = mod1(i + div(nλ, 2), nλ)

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

@inline function _limited_slope(center, backward, forward)
    raw = 0.5 * (forward - backward)
    lower = min(backward, center, forward)
    upper = max(backward, center, forward)
    return copysign(min(abs(raw), 2.0 * (center - lower), 2.0 * (upper - center)), raw)
end

function _face_velocities!(workspace, u, v)
    nλ, nθ, nd = size(u)
    uf, vf = workspace.u_face, workspace.v_face
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        uf[i, j, k] = 0.5 * (u[_west(i, nλ), j, k] + u[i, j, k])
    end
    fill!(vf, 0.0)
    @inbounds for k = 1:nd, j = 2:nθ, i = 1:nλ
        vf[i, j, k] = 0.5 * (v[i, j-1, k] + v[i, j, k])
    end
    return nothing
end

function _horizontal_substeps(workspace, mesh, dt, cfl_limit)
    max_cfl = 0.0
    nλ, nθ, nd = size(workspace.u_face)
    Δλ = 2π / nλ
    uf, vf = workspace.u_face, workspace.v_face
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        ie = _east(i, nλ)
        zonal_length = mesh.radius * (mesh.θe[j+1] - mesh.θe[j])
        west_volume = uf[i, j, k] * zonal_length
        east_volume = uf[ie, j, k] * zonal_length
        south_volume = vf[i, j, k] * mesh.radius * cos(mesh.θe[j]) * Δλ
        north_volume = vf[i, j+1, k] * mesh.radius * cos(mesh.θe[j+1]) * Δλ
        incoming = max(west_volume, 0.0) + max(-east_volume, 0.0) +
                   max(south_volume, 0.0) + max(-north_volume, 0.0)
        area = mesh.radius^2 * Δλ * (sin(mesh.θe[j+1]) - sin(mesh.θe[j]))
        max_cfl = max(max_cfl, dt * incoming / area)
    end
    return max(1, ceil(Int, max_cfl / cfl_limit))
end

function _vertical_substeps(M, Δp, dt, cfl_limit)
    nλ, nθ, nd = size(Δp)
    max_cfl = 0.0
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        incoming = max(M[i, j, k], 0.0) + max(-M[i, j, k+1], 0.0)
        max_cfl = max(max_cfl, dt * incoming / Δp[i, j, k])
    end
    return max(1, ceil(Int, max_cfl / cfl_limit))
end

"""
Advance a scalar with spherical multidimensional Van Leer transport followed by
nonuniform monotone PPM vertical transport. Winds, mass fluxes, and layer
thicknesses are held fixed while the two operators subcycle independently.
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
    all(isfinite, tracer_in) && minimum(tracer_in) >= 0.0 ||
        throw(DomainError(minimum(tracer_in), "input tracer must be finite and non-negative"))
    all(isfinite, Δp) && minimum(Δp) > 0.0 ||
        throw(DomainError(minimum(Δp), "layer pressure thickness must be finite and positive"))

    _face_velocities!(workspace, u, v)
    copyto!(workspace.q_work, tracer_in)

    nh = _horizontal_substeps(workspace, mesh, dt, cfl_limit)
    dt_h = dt / nh
    for _ = 1:nh
        Horizontal_Tracer_Advection!(workspace, mesh, workspace.q_stage, workspace.q_work, u, v, dt_h)
        _remove_roundoff_undershoots!(workspace.q_stage, "horizontal")
        workspace.q_work, workspace.q_stage = workspace.q_stage, workspace.q_work
    end

    nv = _vertical_substeps(M, Δp, dt, cfl_limit)
    dt_v = dt / nv
    for _ = 1:nv
        Vertical_Tracer_Advection!(workspace, workspace.q_stage, workspace.q_work, Δp, M, dt_v)
        _remove_roundoff_undershoots!(workspace.q_stage, "vertical")
        workspace.q_work, workspace.q_stage = workspace.q_stage, workspace.q_work
    end

    copyto!(tracer_out, workspace.q_work)
    all(isfinite, tracer_out) || error("grid tracer transport produced non-finite values")
    minimum(tracer_out) >= 0.0 || error("grid tracer roundoff cleanup failed")
    return (horizontal_substeps = nh, vertical_substeps = nv)
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
        c = 0.5dt * u[i, j, k] / (mesh.radius * mesh.cosθ[j] * Δλ)
        qx[i, j, k] = c >= 0.0 ?
            q[i, j, k] + c * (q[_west(i, nλ), j, k] - q[i, j, k]) :
            q[i, j, k] - c * (q[_east(i, nλ), j, k] - q[i, j, k])

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

"""One CFL-limited multidimensional spherical Van Leer substep."""
function Horizontal_Tracer_Advection!(workspace, mesh, q_out, q, u, v, dt)
    nλ, nθ, nd = size(q)
    Δλ = 2π / nλ
    _transverse_states!(workspace, mesh, q, u, v, dt)
    _longitude_slopes!(workspace.slope_lambda, workspace.q_half_phi)
    _latitude_slopes!(workspace.slope_phi, workspace.q_half_lambda, mesh)

    uf, vf = workspace.u_face, workspace.v_face
    fqλ, flλ, fvλ = workspace.tracer_flux_lambda,
    workspace.low_flux_lambda,
    workspace.volume_flux_lambda
    fqφ, flφ, fvφ = workspace.tracer_flux_phi,
    workspace.low_flux_phi,
    workspace.volume_flux_phi

    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        im = _west(i, nλ)
        if uf[i, j, k] >= 0.0
            donor = im
            c = uf[i, j, k] * dt / (mesh.radius * mesh.cosθ[j] * Δλ)
            qface = workspace.q_half_phi[donor, j, k] + 0.5 * workspace.slope_lambda[donor, j, k] * (1 - c)
        else
            donor = i
            c = -uf[i, j, k] * dt / (mesh.radius * mesh.cosθ[j] * Δλ)
            qface = workspace.q_half_phi[donor, j, k] - 0.5 * workspace.slope_lambda[donor, j, k] * (1 - c)
        end
        face_length = mesh.radius * (mesh.θe[j+1] - mesh.θe[j])
        fvλ[i, j, k] = uf[i, j, k] * face_length
        fqλ[i, j, k] = fvλ[i, j, k] * qface
        flλ[i, j, k] = fvλ[i, j, k] * q[donor, j, k]
    end

    fill!(fqφ, 0.0)
    fill!(flφ, 0.0)
    fill!(fvφ, 0.0)
    @inbounds for k = 1:nd, h = 2:nθ, i = 1:nλ
        if vf[i, h, k] >= 0.0
            donor = h - 1
            c = vf[i, h, k] * dt / (mesh.radius * (mesh.θe[donor+1] - mesh.θe[donor]))
            qface = workspace.q_half_lambda[i, donor, k] + 0.5 * workspace.slope_phi[i, donor, k] * (1 - c)
        else
            donor = h
            c = -vf[i, h, k] * dt / (mesh.radius * (mesh.θe[donor+1] - mesh.θe[donor]))
            qface = workspace.q_half_lambda[i, donor, k] - 0.5 * workspace.slope_phi[i, donor, k] * (1 - c)
        end
        face_length = mesh.radius * cos(mesh.θe[h]) * Δλ
        fvφ[i, h, k] = vf[i, h, k] * face_length
        fqφ[i, h, k] = fvφ[i, h, k] * qface
        flφ[i, h, k] = fvφ[i, h, k] * q[i, donor, k]
    end

    # Flux-corrected transport: the donor-cell update is positive under the
    # combined incoming-CFL bound. Limit only the high-order antidiffusive
    # correction, using one factor per losing cell and one shared face flux.
    ratio = workspace.positivity_ratio
    tolerance = TRACER_ROUNDOFF_FACTOR * eps(max(maximum(abs, q), 1.0))
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        ie = _east(i, nλ)
        area = mesh.radius^2 * Δλ * (sin(mesh.θe[j+1]) - sin(mesh.θe[j]))
        low_div = flλ[ie, j, k] - flλ[i, j, k] + flφ[i, j+1, k] - flφ[i, j, k]
        volume_div = fvλ[ie, j, k] - fvλ[i, j, k] + fvφ[i, j+1, k] - fvφ[i, j, k]
        q_low = q[i, j, k] + dt * (-low_div + q[i, j, k] * volume_div) / area
        q_low >= -tolerance || error("horizontal donor-cell update violated positivity: minimum=$q_low")

        correction_e = fqλ[ie, j, k] - flλ[ie, j, k]
        correction_w = fqλ[i, j, k] - flλ[i, j, k]
        correction_n = fqφ[i, j+1, k] - flφ[i, j+1, k]
        correction_s = fqφ[i, j, k] - flφ[i, j, k]
        loss = max(correction_e, 0.0) + max(-correction_w, 0.0) +
               max(correction_n, 0.0) + max(-correction_s, 0.0)
        ratio[i, j, k] = loss > 0.0 ? min(1.0, max(q_low, 0.0) * area / (dt * loss)) : 1.0
    end

    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        correction = fqλ[i, j, k] - flλ[i, j, k]
        losing_i = correction >= 0.0 ? _west(i, nλ) : i
        fqλ[i, j, k] = flλ[i, j, k] + ratio[losing_i, j, k] * correction
    end
    @inbounds for k = 1:nd, h = 2:nθ, i = 1:nλ
        correction = fqφ[i, h, k] - flφ[i, h, k]
        losing_j = correction >= 0.0 ? h - 1 : h
        fqφ[i, h, k] = flφ[i, h, k] + ratio[i, losing_j, k] * correction
    end

    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        ie = _east(i, nλ)
        area = mesh.radius^2 * Δλ * (sin(mesh.θe[j+1]) - sin(mesh.θe[j]))
        tracer_div = fqλ[ie, j, k] - fqλ[i, j, k] + fqφ[i, j+1, k] - fqφ[i, j, k]
        volume_div = fvλ[ie, j, k] - fvλ[i, j, k] + fvφ[i, j+1, k] - fvφ[i, j, k]
        tendency = (-tracer_div + q[i, j, k] * volume_div) / area
        workspace.tendency[i, j, k] = tendency
        q_out[i, j, k] = q[i, j, k] + dt * tendency
    end
    return nothing
end

function _ppm_reconstruction!(workspace, q, dz)
    nλ, nθ, nd = size(q)
    s, ql, qr, q6 = workspace.ppm_slope, workspace.q_left, workspace.q_right, workspace.q6
    fill!(s, 0.0)
    @inbounds for j = 1:nθ, i = 1:nλ, k = 2:nd-1
        gm = (q[i, j, k] - q[i, j, k-1]) / (dz[i, j, k] + dz[i, j, k-1])
        gp = (q[i, j, k+1] - q[i, j, k]) / (dz[i, j, k+1] + dz[i, j, k])
        raw = (gp * (2dz[i, j, k-1] + dz[i, j, k]) +
               gm * (2dz[i, j, k+1] + dz[i, j, k])) * dz[i, j, k] /
              (dz[i, j, k-1] + dz[i, j, k] + dz[i, j, k+1])
        lower = min(q[i, j, k-1], q[i, j, k], q[i, j, k+1])
        upper = max(q[i, j, k-1], q[i, j, k], q[i, j, k+1])
        s[i, j, k] = copysign(min(abs(raw), 2(q[i, j, k] - lower), 2(upper - q[i, j, k])), raw)
    end
    @. ql = q - 0.5s
    @. qr = q + 0.5s

    @inbounds for j = 1:nθ, i = 1:nλ, h = 3:nd-1
        km, kp = h - 1, h
        hmm, hm, hp, hpp = dz[i, j, km-1], dz[i, j, km], dz[i, j, kp], dz[i, j, kp+1]
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

    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
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
    return nothing
end

"""One CFL-limited monotone nonuniform-PPM vertical substep."""
function Vertical_Tracer_Advection!(workspace, q_out, q, Δp, M, dt)
    nλ, nθ, nd = size(q)
    _ppm_reconstruction!(workspace, q, Δp)
    flux, low_flux = workspace.vertical_flux, workspace.low_vertical_flux
    fill!(flux, 0.0)
    fill!(low_flux, 0.0)
    @inbounds for h = 2:nd, j = 1:nθ, i = 1:nλ
        if M[i, j, h] >= 0.0
            donor = h - 1
            c = M[i, j, h] * dt / Δp[i, j, donor]
            swept = workspace.q_right[i, j, donor] - 0.5c * (
                workspace.q_right[i, j, donor] - workspace.q_left[i, j, donor] -
                (1 - 2c / 3) * workspace.q6[i, j, donor]
            )
        else
            donor = h
            c = -M[i, j, h] * dt / Δp[i, j, donor]
            swept = workspace.q_left[i, j, donor] + 0.5c * (
                workspace.q_right[i, j, donor] - workspace.q_left[i, j, donor] +
                (1 - 2c / 3) * workspace.q6[i, j, donor]
            )
        end
        flux[i, j, h] = M[i, j, h] * swept
        low_flux[i, j, h] = M[i, j, h] * q[i, j, donor]
    end

    ratio = workspace.positivity_ratio
    tolerance = TRACER_ROUNDOFF_FACTOR * eps(max(maximum(abs, q), 1.0))
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        low_tendency = -(low_flux[i, j, k+1] - low_flux[i, j, k] -
                         q[i, j, k] * (M[i, j, k+1] - M[i, j, k])) / Δp[i, j, k]
        q_low = q[i, j, k] + dt * low_tendency
        q_low >= -tolerance || error("vertical donor-cell update violated positivity: minimum=$q_low")
        correction_bottom = flux[i, j, k+1] - low_flux[i, j, k+1]
        correction_top = flux[i, j, k] - low_flux[i, j, k]
        loss = max(correction_bottom, 0.0) + max(-correction_top, 0.0)
        ratio[i, j, k] = loss > 0.0 ?
            min(1.0, max(q_low, 0.0) * Δp[i, j, k] / (dt * loss)) : 1.0
    end
    @inbounds for h = 2:nd, j = 1:nθ, i = 1:nλ
        correction = flux[i, j, h] - low_flux[i, j, h]
        losing_k = correction >= 0.0 ? h - 1 : h
        flux[i, j, h] = low_flux[i, j, h] + ratio[i, j, losing_k] * correction
    end
    @inbounds for k = 1:nd, j = 1:nθ, i = 1:nλ
        tendency = -(flux[i, j, k+1] - flux[i, j, k] -
                     q[i, j, k] * (M[i, j, k+1] - M[i, j, k])) / Δp[i, j, k]
        workspace.tendency[i, j, k] = tendency
        q_out[i, j, k] = q[i, j, k] + dt * tendency
    end
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
