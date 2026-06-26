module Vertical_Interpolation_Module

using Base.Threads
using ..Atmo_Data_Module

export Compute_Pressure_Grid!, Interpolate_Field!



function Compute_Pressure_Grid!(
    p_3d::AbstractArray{Float64,3},
    ak::Vector{Float64},
    bk::Vector{Float64},
    ps::AbstractArray{Float64,2},
)
    nλ, nθ, nd = size(p_3d)

    # Parallelize over latitude (outer loop) to minimize overhead
    @threads for j = 1:nθ
        for i = 1:nλ
            @inbounds ps_val = ps[i, j]
            @inbounds for k = 1:nd
                p_3d[i, j, k] = ak[k] + bk[k] * ps_val
            end
        end
    end
end

using Base.Threads



function Interpolate_Field!(
    out_3d::AbstractArray{Float64,3},
    in_3d::AbstractArray{Float64,3},
    p_3d::AbstractArray{Float64,3},
    ps_2d::AbstractArray{Float64,2},
    log_targets::Vector{Float64},
    var_name::Symbol,
    atmo_data::Atmo_Data,
    t_3d::Union{AbstractArray{Float64,3},Nothing} = nothing,
)
    nλ, nθ, n_plev = size(out_3d)
    nd = size(in_3d, 3)

    # 1. Extract Constants
    R_d = atmo_data.rdgas
    ALPHA = atmo_data.alpha

    # 2. Determine Strategy
    strategy = :constant
    if var_name in (:t, :t_eq)
        strategy = :lapse_rate
    elseif var_name == :z
        strategy = :hydrostatic
    end

    @threads for j = 1:nθ
        for i = 1:nλ
            # Lightweight views
            col_in = view(in_3d, i, j, :)
            col_p = view(p_3d, i, j, :)

            # True Surface Pressure for this column (from prognostic variable)
            ps_true = ps_2d[i, j]

            # For Hydrostatic strategy, we need the T column
            col_t =
                (strategy == :hydrostatic && !isnothing(t_3d)) ? view(t_3d, i, j, :) :
                nothing

            # Check direction (Standard: p[1] is Top, p[end] is Bottom)
            is_down = col_p[end] > col_p[1]

            for k_tgt = 1:n_plev
                lt = log_targets[k_tgt]
                p_target = exp(lt)

                # --- EXTRAPOLATION (Underground) ---
                # Strictly check against True Surface Pressure
                if p_target > ps_true
                    # Variables at the Lowest Model Level (e.g. sigma = 0.99)
                    val_lowest_level = col_in[end]
                    p_lowest_level = col_p[end]

                    if strategy == :constant
                        out_3d[i, j, k_tgt] = val_lowest_level

                    elseif strategy == :lapse_rate
                        # Extrapolate T from the lowest model level down to target
                        out_3d[i, j, k_tgt] =
                            val_lowest_level * (p_target / p_lowest_level)^ALPHA

                    elseif strategy == :hydrostatic
                        # Geopotential Extrapolation
                        if isnothing(col_t)
                            out_3d[i, j, k_tgt] = NaN
                        else
                            T_lowest = col_t[end]
                            T_tgt = T_lowest * (p_target / p_lowest_level)^ALPHA
                            T_mean = 0.5 * (T_lowest + T_tgt)

                            # Hydrostatic Integral from lowest model level to target
                            # Phi_tgt = Phi_lowest - R*T_mean*ln(P_tgt/P_lowest)
                            out_3d[i, j, k_tgt] =
                                val_lowest_level -
                                (R_d * T_mean) * (lt - log(p_lowest_level))
                        end
                    else
                        out_3d[i, j, k_tgt] = val_lowest_level
                    end
                    continue
                end

                # --- INTERPOLATION (In Air) ---
                min_lp = log(is_down ? col_p[1] : col_p[end])

                if lt < min_lp # Above Top
                    out_3d[i, j, k_tgt] = NaN
                    continue
                end

                k_found = -1
                for k = 1:nd-1
                    p_k = max(col_p[k], 1e-10)
                    p_kp1 = max(col_p[k+1], 1e-10)
                    lp1 = log(p_k)
                    lp2 = log(p_kp1)
                    if (lp1 <= lt <= lp2) || (lp2 <= lt <= lp1)
                        k_found = k
                        break
                    end
                end

                if k_found != -1
                    p_k = max(col_p[k_found], 1e-10)
                    p_kp1 = max(col_p[k_found+1], 1e-10)
                    lp1 = log(p_k)
                    lp2 = log(p_kp1)
                    v1 = col_in[k_found]
                    v2 = col_in[k_found+1]
                    w = (lt - lp1) / (lp2 - lp1)
                    out_3d[i, j, k_tgt] = v1 + w * (v2 - v1)
                else
                    # Fallback for the thin layer between P_lowest_level and Ps_true
                    # (Where p_target < Ps_true but > P_lowest_level)
                    out_3d[i, j, k_tgt] = col_in[end]
                end
            end
        end
    end
end

end
