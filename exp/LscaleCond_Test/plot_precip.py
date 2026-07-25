"""
Plot time-mean precipitation map from LscaleCond_Test output.
Checks that precipitation is confined to the tropics (|lat| <= 30 deg).

Native output is written as time-stamped chunks (output_t0.nc, output_t<N>.nc, ...)
rather than a single output.nc — see Output_Manager's chunked NC output strategy.
Pressure-level files (output_t*_plev.nc) are excluded; this script only needs
the native-grid "pr" field. netCDF4.MFDataset can't be used to join the chunks:
NCDatasets.jl writes full NETCDF4 (HDF5) files, and MFDataset only supports the
NETCDF3_*/NETCDF4_CLASSIC formats — xarray.open_mfdataset has no such restriction.
"""
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import xarray as xr

output_pattern = "output_t*.nc"
plot_file      = "precip_map.png"

chunk_files = sorted(f for f in glob.glob(output_pattern) if "_plev" not in f)
if not chunk_files:
    raise FileNotFoundError(f"No files matching '{output_pattern}' found in current directory")

print(f"Found {len(chunk_files)} output chunk(s): {chunk_files}")

with xr.open_mfdataset(
    chunk_files, combine="by_coords",
    data_vars="minimal", coords="minimal", compat="override"
) as ds:
    lon = ds["lon"].values   # [nλ]
    lat = ds["lat"].values   # [nθ]

    precip_var = ds["pr"]
    precip = precip_var.values   # [time, lat, lon]

    print("precip shape:", precip.shape)
    print("precip units:", precip_var.attrs.get("units", "unknown"))

# Time-mean
precip_mean = precip.mean(axis=0)   # [lat, lon] or [lon, lat]

# Ensure shape is [lat, lon]
if precip_mean.shape[0] == len(lon) and precip_mean.shape[1] == len(lat):
    precip_mean = precip_mean.T     # transpose to [lat, lon]

print(f"lat range : {lat.min():.1f} to {lat.max():.1f}")
print(f"precip min/max: {precip_mean.min():.4e} / {precip_mean.max():.4e}")

# ---- Plot ----
fig, ax = plt.subplots(figsize=(12, 5))

LON, LAT = np.meshgrid(lon, lat)
cf = ax.contourf(LON, LAT, precip_mean, levels=20, cmap="Blues")
cb = plt.colorbar(cf, ax=ax, label="Precipitation (time-mean)")

# 30S / 30N reference lines
ax.axhline( 30, color="red",  lw=1.5, ls="--", label="±30°")
ax.axhline(-30, color="red",  lw=1.5, ls="--")

ax.set_xlabel("Longitude (°E)")
ax.set_ylabel("Latitude (°N)")
ax.set_title("LscaleCond_Test — Time-mean precipitation\n"
             "(L=0.5 in tropics, L=0 elsewhere)")
ax.legend(loc="upper right")
ax.set_yticks(np.arange(-90, 91, 30))
ax.set_xticks(np.arange(0, 361, 60))

plt.tight_layout()
plt.savefig(plot_file, dpi=150)
print(f"Saved: {plot_file}")
plt.show()
