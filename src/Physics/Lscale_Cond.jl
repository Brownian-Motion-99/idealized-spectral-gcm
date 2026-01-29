using Base.Threads
using ...Vert_Coordinate_Module
using ...Atmo_Data_Module

function Lscale_Cond!(
    vert_coord::Vert_Coordinate,
    atmo_data::Atmo_Data, 
    grid_q::AbstractArray{Float64, 3}, grid_δq::AbstractArray{Float64, 3}, grid_liquid_water_content::Array{Float64, 3}, grid_precip::Array{Float64, 3},
    grid_t::Array{Float64, 3}, grid_δt::Array{Float64, 3}, 
    grid_p_full::Array{Float64, 3}, grid_ps::Array{Float64, 3}, 
    Δt
)

    Δak, Δbk   = vert_coord.Δak, vert_coord.Δbk
    nλ, nθ, nd = atmo_data.nλ, atmo_data.nθ, atmo_data.nd

    cp = atmo_data.cp_air
    Lv = atmo_data.Lv
    Rv = atmo_data.rvgas
    L  = atmo_data.L

    const_es     = 611.12
    const_q1     = 0.622
    const_q2     = 0.378
    Lv_Rv        = Lv / Rv
    inv_273      = 1.0 / 273.15
    Lv_cp        = Lv / cp
    heating_rate = L * Lv_cp

    @threads for k = 1:nd
        for j = 1:nθ
            for i = 1:nλ

                # Load data
                t_val = grid_t[i, j, k]
                p_val = grid_p_full[i, j, k]
                q_val = grid_q[i, j, k]

                # Saturated specific humidity
                es     = const_es * exp(Lv_Rv * (inv_273 - 1.0 / t_val))
                qs     = (const_q1 * es) / (p_val - const_q2 * es)
                dqs_dt = Lv * qs / (Rv * t_val^2)

                # Heating
                δq               = (max(q_val, qs) - qs) / (1.0 + Lv_cp * dqs_dt) / (2.0 * Float64(Δt))
                grid_δq[i, j, k] = δq
                grid_δt[i, j, k] = δq * heating_rate

                # Liquid water content
                grid_liquid_water_content[i, j, k] = δq

                # Pseudo-adiabatic precipitation
                Δp                    = Δak[k] + Δbk[k] * grid_ps[i, j, 1]
                grid_precip[i, j, 1] += δq * Δp

            end
        end
    end

end
