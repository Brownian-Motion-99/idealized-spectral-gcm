"""
Plot time-mean precipitation map from LscaleCond_Test output.
Checks that precipitation is confined to the tropics (|lat| <= 30 deg).
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import netCDF4 as nc

output_file = "exp/LscaleCond_Test/output.nc"
plot_file   = "exp/LscaleCond_Test/precip_map.png"

with nc.Dataset(output_file) as ds:
    # Axes
    lon = ds.variables["lon"][:]   # [nλ]
    lat = ds.variables["lat"][:]   # [nθ]

    # precipitation [time, lat, lon] or [time, lon, lat] — detect
    precip_var = ds.variables["pr"]
    precip = precip_var[:]         # load all timesteps

    print("precip shape:", precip.shape)
    print("precip units:", getattr(precip_var, "units", "unknown"))

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
