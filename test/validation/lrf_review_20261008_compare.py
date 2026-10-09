"""Summarize the paired BM-control/LRF pilot's saved interval means.

Usage: python lrf_review_20261008_compare.py CONTROL_NC LRF_NC OUTPUT_JSON
Column heating uses nominal layer pressures and Gaussian area weights. It is
an approximate heating diagnostic, not an exact instantaneous energy ledger.
"""
import json
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def compare(control_path, treatment_path):
    with Dataset(control_path) as control, Dataset(treatment_path) as treatment:
        for name in ("time", "lat", "lon", "pfull", "phalf"):
            np.testing.assert_array_equal(control[name][:], treatment[name][:])
        assert treatment["lrf_dta_dt"].dimensions == ("time", "pfull", "lat", "lon")
        heat = np.asarray(treatment["lrf_dta_dt"][:])
        assert np.isfinite(heat).all()
        assert np.all(control["lrf_dta_dt"][:] == 0)
        assert np.any(heat != 0)
        nodes, weights = np.polynomial.legendre.leggauss(len(treatment["lat"]))
        gaussian_latitude = np.rad2deg(np.arcsin(nodes))
        latitude = np.asarray(treatment["lat"][:])
        indices = np.abs(latitude[:, None] - gaussian_latitude).argmin(axis=1)
        np.testing.assert_allclose(latitude, gaussian_latitude[indices], atol=1e-10, rtol=0)
        assert len(np.unique(indices)) == len(indices)
        area_weights = weights[indices] / 2
        dp = np.diff(treatment["phalf"][:]) * 100
        cp = float(treatment.dry_air_gas_constant) / (2 / 7)
        gravity = float(treatment.gravity)
        column_heating = np.einsum("tklm,k,l->t", heat, dp, area_weights)
        column_heating *= cp / (gravity * heat.shape[-1])
        mean_profile = heat.mean(axis=(0, 3)) * 86400
        profile_rms = np.sqrt(np.einsum("kl,k,l->", mean_profile**2, dp, area_weights) / dp.sum())
        differences = {}
        for name in ("ta", "hus", "ua", "va", "ps"):
            difference = np.asarray(treatment[name][-1]) - np.asarray(control[name][-1])
            assert np.isfinite(difference).all()
            differences[name] = {
                "maximum_absolute": float(np.abs(difference).max()),
                "rms": float(np.sqrt(np.mean(difference**2))),
            }
        return {
            "control": str(Path(control_path).resolve()),
            "treatment": str(Path(treatment_path).resolve()),
            "times_day": treatment["time"][:].tolist(),
            "lrf_grid_rms_K_day": (np.sqrt(np.mean(heat**2, axis=(1, 2, 3))) * 86400).tolist(),
            "nominal_global_column_heating_W_m2": column_heating.tolist(),
            "one_day_time_zonal_mean_profile_rms_K_day": float(profile_rms),
            "last_interval_mean_differences": differences,
            "specific_heat_J_kg_K": cp,
            "gravity_m_s2": gravity,
            "weighting": "Gaussian horizontal area; nominal pressure layer mass",
            "scope": "saved interval means; first day only; not equilibrium climate validation",
        }


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    result = compare(sys.argv[1], sys.argv[2])
    Path(sys.argv[3]).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
