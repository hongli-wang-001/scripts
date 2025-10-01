import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import cartopy.mpl.ticker as cticker
import os

# === Configuration ===
n_ens = 30
ens_dir = "data/ens"
ens_template = "mem{:03d}/mpasout.2024-05-27_00.00.00.nc"
invariant_file = "invariant.nc"
output_file = "ens_mean.nc"
level_index = 23  # used for plotting only
variables_to_process = ['surface_pressure', 'uReconstructZonal', 'uReconstructMeridional', 'pressure_p', 'theta']

# === Load Static Grid ===
with Dataset(invariant_file, 'r') as nc:
    latCell = np.degrees(nc.variables['latCell'][:])
    lonCell = np.degrees(nc.variables['lonCell'][:])
    terrain = nc.variables['ter'][:]

lat_min, lat_max = latCell.min() - 2, latCell.max() + 2
lon_min, lon_max = lonCell.min() - 2, lonCell.max() + 2

# === Determine Dimensions from First Ensemble Member ===
sample_file = os.path.join(ens_dir, ens_template.format(1))
with Dataset(sample_file, 'r') as nc:
    nCells = len(nc.dimensions['nCells'])
    nVertLevels = len(nc.dimensions['nVertLevels'])

# === Read and Collect Ensemble Data ===
ens_data = {}
for var in variables_to_process:
    shape = (n_ens, nCells, nVertLevels) if var != 'surface_pressure' else (n_ens, nCells)
    ens_data[var] = np.zeros(shape, dtype=np.float32)

for i in range(n_ens):
    filepath = os.path.join(ens_dir, ens_template.format(i + 1))
    with Dataset(filepath, 'r') as nc:
        for var in variables_to_process:
            data = nc.variables[var]
            if data.ndim == 3:
                ens_data[var][i, :, :] = data[0, :, :]
            else:
                ens_data[var][i, :] = data[0, :]

# === Compute Ensemble Mean and Std Dev ===
ens_stats = {}
for var in variables_to_process:
    mean = np.mean(ens_data[var], axis=0)
    std = np.std(ens_data[var], axis=0)
    ens_stats[var] = {'mean': mean, 'std': std}

# === Plot Function ===
def plot_var(var_data, title, filename, cmap='viridis'):
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_title(title, fontsize=14)
    ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.LAND, facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN)

    gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = cticker.LongitudeFormatter()
    gl.yformatter = cticker.LatitudeFormatter()

    sc = ax.scatter(lonCell, latCell, c=var_data, cmap=cmap, s=1, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(sc, orientation='vertical', pad=0.02, shrink=0.75)
    cbar.set_label(title)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()
    print(f"✅ Saved: {filename}")

# === Plot Terrain ===
plot_var(terrain, "Terrain Height (m)", "terrain.png", cmap='terrain')

# === Plot Ensemble Mean and Std Dev at level_index ===
for var in variables_to_process:
    if ens_stats[var]['mean'].ndim == 2:
        plot_var(ens_stats[var]['mean'][:, level_index], f"{var} Mean (level={level_index})", f"{var}_mean.png")
        plot_var(ens_stats[var]['std'][:, level_index], f"{var} Std Dev (level={level_index})", f"{var}_std.png")
    else:
        plot_var(ens_stats[var]['mean'], f"{var} Mean", f"{var}_mean.png")
        plot_var(ens_stats[var]['std'], f"{var} Std Dev", f"{var}_std.png")

# === Write Ensemble Mean and Std to NetCDF ===
with Dataset(output_file, 'w') as nc_out:
    # Dimensions
    nc_out.createDimension('nCells', nCells)
    nc_out.createDimension('nVertLevels', nVertLevels)

    # Coordinates
    lat_var = nc_out.createVariable('latCell', 'f4', ('nCells',))
    lon_var = nc_out.createVariable('lonCell', 'f4', ('nCells',))
    lat_var[:] = latCell
    lon_var[:] = lonCell
    lat_var.units = 'degrees_north'
    lon_var.units = 'degrees_east'

    # Write each variable
    for var in variables_to_process:
        if var == 'surface_pressure':
            shape = ('nCells',)
            nc_mean = nc_out.createVariable(f"{var}_mean", 'f4', shape, zlib=True)
            nc_std = nc_out.createVariable(f"{var}_std", 'f4', shape, zlib=True)
            nc_mean[:] = ens_stats[var]['mean']
            nc_std[:] = ens_stats[var]['std']
        else:
            shape = ('nCells', 'nVertLevels')
            nc_mean = nc_out.createVariable(f"{var}_mean", 'f4', shape, zlib=True)
            nc_std = nc_out.createVariable(f"{var}_std", 'f4', shape, zlib=True)
            nc_mean[:, :] = ens_stats[var]['mean']
            nc_std[:, :] = ens_stats[var]['std']

        nc_mean.long_name = f"{var} ensemble mean"
        nc_std.long_name = f"{var} ensemble std dev"

print(f"✅ NetCDF ensemble mean/std written to: {output_file}")

