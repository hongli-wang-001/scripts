import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
from matplotlib.colors import TwoSlopeNorm

# Enlarge all font sizes globally
plt.rcParams.update({
    'font.size': 14,            # Default text size
    'axes.titlesize': 16,       # Axes title
    'axes.labelsize': 14,       # Axes labels
    'xtick.labelsize': 12,      # Tick labels
    'ytick.labelsize': 12,
    'legend.fontsize': 12,      # Legend
    'figure.titlesize': 16      # Figure title (if used)
})


# ============================
# USER CONFIGURATION
# ============================
nc_path = 'combined_ioda_v29.nc'     # Path to input NetCDF file
channel_index = 59 #141 #0 #94                # Radiance channel (0-based index)
domain_extent = [-90, 90, 0, 89]     # [lon_min, lon_max, lat_min, lat_max]
cmap = "jet"                         # Colormap (e.g., "inferno", "viridis", "RdBu_r")
output_figure = f"radiance_chan{channel_index + 1}.png"
# ============================

# --- Open NetCDF ---
ds = Dataset(nc_path, mode='r')
meta = ds.groups['MetaData']
obs = ds.groups['ObsValue']

# --- Helper function for masking fill values ---
def mask_fill(var):
    fill = getattr(var, '_FillValue', None)
    data = var[:]
    return np.ma.masked_where(data == fill, data) if fill is not None else data

# --- Read and mask variables ---
lat = mask_fill(meta.variables['latitude'])
lon = mask_fill(meta.variables['longitude'])
zenith = mask_fill(meta.variables['sensorZenithAngle'])
azimuth = mask_fill(meta.variables['sensorAzimuthAngle'])
cloud = mask_fill(meta.variables['cloudAmount'])

# Handle invalid placeholders (e.g., 65540.0)
invalid_val = 65540.0
mask_invalid = (lat == invalid_val) | (lon == invalid_val)
lat = np.ma.masked_where(mask_invalid, lat)
lon = np.ma.masked_where(mask_invalid, lon)
zenith = np.ma.masked_where(mask_invalid, zenith)

# --- Read and mask radiance for selected channel ---
rec_rad = obs.variables['radiance']
radiance = rec_rad[:, channel_index]

# Mask out fill values and physical bounds
fill_value = getattr(rec_rad, '_FillValue', None)
range_mask = (radiance < 0.0) | (radiance > 20.0)
if fill_value is not None:
    range_mask |= (radiance == fill_value)

total_mask = mask_invalid | range_mask
radiance = np.ma.masked_where(total_mask, radiance)

# ============================
# PLOTTING FUNCTION
# ============================
def plot_radiance_map(lon, lat, radiance, extent, channel_index, cmap, output_filename):
    lon_min, lon_max, lat_min, lat_max = extent

    # Mask for plot region
    extent_mask = (
        (lon >= lon_min) & (lon <= lon_max) &
        (lat >= lat_min) & (lat <= lat_max)
    )
    radiance_plot = np.ma.masked_where(~extent_mask, radiance)

    # Color normalization
    vmin = radiance_plot.min()
    vmax = radiance_plot.max()
    vcenter = 0.5 * (vmin + vmax)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    print(f"\n📊 Radiance stats for Channel {channel_index + 1} in extent:")
    print(f"   Min    = {vmin:.4f}")
    print(f"   Max    = {vmax:.4f}")
    print(f"   Center = {vcenter:.4f}")

    # --- Create plot ---
    fig = plt.figure(figsize=(12, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.coastlines('50m')

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    # Scatter plot
    sc = ax.scatter(
        lon, lat, c=radiance, cmap=cmap, norm=norm,
        s=1, alpha=0.5, transform=ccrs.PlateCarree()
    )

    # Colorbar
    cbar = plt.colorbar(sc, orientation='horizontal', shrink=0.6, pad=0.03)
    cbar.set_label('Radiance (W/m²·sr·m⁻¹)')

    plt.title(f'Reconstructed Radiance – Channel {channel_index + 1}', fontsize=14)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    plt.close()

    print(f"✅ Plot saved as: {output_filename}")

# ============================
# CALL THE FUNCTION
# ============================
plot_radiance_map(
    lon=lon,
    lat=lat,
    radiance=radiance,
    extent=domain_extent,
    channel_index=channel_index,
    cmap=cmap,
    output_filename=output_figure
)

# --- Close dataset ---
ds.close()

