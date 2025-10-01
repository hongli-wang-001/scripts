import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
from scipy.fftpack import dct, idct
import matplotlib.pyplot as plt
from pyproj import Proj, Transformer
import shutil

# === Configuration ===
remove_ens_mean = True  # Set to False to skip ensemble mean removal
level_index = 1
nx, ny = 224, 134
dx, dy = 35000, 35000 # meter
lon0, lat0 = 261, 41  # center in degrees
invariant_file = 'invariant.nc'
mpas_file = 'mpasout.nc'
mpas_mean_file = 'ens_mean.nc'
output_nc_file = 'interpolated_lambert.nc'
mpas_copy_a = 'mpasout_scale_a.nc'
mpas_copy_b = 'mpasout_scale_b.nc'

# Lambert Conformal Projection
proj_lambert = Proj(proj='lcc', lat_1=lat0, lat_2=lat0, lat_0=lat0,
                    lon_0=lon0, x_0=0, y_0=0, ellps='WGS84')
transformer = Transformer.from_proj("epsg:4326", proj_lambert, always_xy=True)
inv_transformer = Transformer.from_proj(proj_lambert, "epsg:4326", always_xy=True)

# === Read MPAS lat/lon ===
with Dataset(invariant_file, 'r') as nc:
    lat_1d = np.degrees(nc.variables['latCell'][:])
    lon_1d = np.degrees(nc.variables['lonCell'][:])
    x_1d, y_1d = transformer.transform(lon_1d, lat_1d)

# === Read MPAS variables ===
with Dataset(mpas_file, 'r') as nc:
    surface_pressure = nc.variables['surface_pressure'][0, :]
    uReconstructZonal = nc.variables['uReconstructZonal'][0, :, :]
    uReconstructMeridional = nc.variables['uReconstructMeridional'][0, :, :]
    pressure_p = nc.variables['pressure_p'][0, :, :]
    theta = nc.variables['theta'][0, :, :]
    qv = nc.variables['qv'][0, :, :]
nVertLevels = theta.shape[1]

# === Read MPAS variables ===
if remove_ens_mean:
    with Dataset(mpas_mean_file, 'r') as ncb:
        surface_pressure_mean = ncb.variables['surface_pressure_mean'][:]
        uReconstructZonal_mean = ncb.variables['uReconstructZonal_mean'][:, :]
        uReconstructMeridional_mean = ncb.variables['uReconstructMeridional_mean'][:, :]
        pressure_p_mean = ncb.variables['pressure_p_mean'][:, :]
        theta_mean = ncb.variables['theta_mean'][:, :]
        qv_mean = ncb.variables['qv_mean'][:, :]
    surface_pressure = surface_pressure - surface_pressure_mean
    uReconstructZonal = uReconstructZonal - uReconstructZonal_mean
    uReconstructMeridional = uReconstructMeridional - uReconstructMeridional_mean
    pressure_p = pressure_p - pressure_p_mean
    theta = theta - theta_mean
    qv = qv - qv_mean

# === Create Lambert Grid ===
x_center, y_center = transformer.transform(lon0, lat0)
x_min = x_center - (nx // 2) * dx
y_min = y_center - (ny // 2) * dy
x_grid = x_min + dx * np.arange(nx)
y_grid = y_min + dy * np.arange(ny)
x_2d, y_2d = np.meshgrid(x_grid, y_grid)

dct_filter=True
#keep_fraction=0.1 # Smaller = more smoothing 0.1 means keep 10% of frequencies
min_scale_km=300
def dct2(a):
    return dct(dct(a.T, norm='ortho').T, norm='ortho')

def idct2(a):
    return idct(idct(a.T, norm='ortho').T, norm='ortho')

def apply_dct_filter(grid, keep_fraction=0.1):
    # Perform 2D DCT
    dct_grid = dct2(grid)

    # Zero out high-frequency components
    nx, ny = dct_grid.shape
    fx, fy = int(keep_fraction * nx), int(keep_fraction * ny)
    dct_grid[fx:, :] = 0
    dct_grid[:, fy:] = 0

    # Reconstruct filtered grid
    return idct2(dct_grid)

def apply_dct_filter_at_scale_a(grid, dx_km, dy_km, min_scale_km):
    nx, ny = grid.shape

    # Compute frequency cutoff index
    fx = int((nx * dx_km) / min_scale_km)
    fy = int((ny * dy_km) / min_scale_km)
    fx = min(fx, fy)
    fy = fx 

    # 2D DCT
    dct_grid = dct2(grid)

    # Zero out small-scale (high-frequency) features
    dct_grid[fx:, :] = 0
    dct_grid[:, fy:] = 0

    # Inverse 2D DCT
    return idct2(dct_grid)

def apply_dct_filter_at_scale(grid, dx_km, dy_km, min_scale_km):
    nx, ny = grid.shape

    # Compute frequency cutoff index
    # DENIS etal DCT paper, MWR, EQ.7
    # For simplicity, here take the square area by nx*dx_km X ny*dy_km.
    fx = int((2 * nx * dx_km) / min_scale_km)
    fy = int((2 * ny * dy_km) / min_scale_km)
    fx = min(fx, fy)
    fy = fx

    # 2D DCT
    dct_grid = dct2(grid)

    # Zero out small-scale (high-frequency) features
    dct_grid[fx:, :] = 0
    dct_grid[:, fy:] = 0

    # Inverse 2D DCT
    return idct2(dct_grid)


# === Interpolation Functions == 
#Uses Delaunay triangulation to build triangles between scattered input points
#For each point on the target grid, it finds the triangle it falls into
#Interpolates linearly between the triangle’s vertices

def interp_irregular_to_regular_a(var_1d):
    points = np.column_stack((x_1d, y_1d))
    return griddata(points, var_1d, (x_2d, y_2d), method='linear')

def interp_irregular_to_regular(var_1d):
    points = np.column_stack((x_1d, y_1d))

    # Step 1: Linear interpolation (may contain NaNs)
    grid_linear = griddata(points, var_1d, (x_2d, y_2d), method='linear')

    # Step 2: Nearest interpolation (for backup)
    grid_nearest = griddata(points, var_1d, (x_2d, y_2d), method='nearest')

    # Step 3: Fill NaNs with nearest
    filled_grid = np.where(np.isnan(grid_linear), grid_nearest, grid_linear)

    # Step 4: Optional DCT smoothing
    if dct_filter:
        #filled_grid = apply_dct_filter(filled_grid, keep_fraction=keep_fraction)
        filled_grid = apply_dct_filter_at_scale(filled_grid, dx_km=35, dy_km=35, min_scale_km=min_scale_km)

    return filled_grid

def interp_3d_to_regular(var_3d):
    result = np.full((nVertLevels, ny, nx), np.nan, dtype=np.float32)
    for k in range(nVertLevels):
        result[k, :, :] = interp_irregular_to_regular(var_3d[:, k])
    return result

def interp_regular_to_irregular(grid_2d):
    return griddata((x_2d.flatten(), y_2d.flatten()), grid_2d.flatten(), (x_1d, y_1d), method='linear')


def interp_3d_to_irregular(grid_3d):
    return np.vstack([interp_regular_to_irregular(grid_3d[k]) for k in range(nVertLevels)]).T

def interp_3d_back(grid_3d):
    return np.vstack([interp_regular_to_irregular(grid_3d[k]) for k in range(nVertLevels)]).T

def update_mpas_file(output_file, surface_pressure_val, u_zonal_val, u_merid_val, pressure_val, theta_val, qv_val):
    shutil.copyfile(mpas_file, output_file)
    with Dataset(output_file, 'r+') as nc:
        nc.variables['surface_pressure'][0, :] = surface_pressure_val
        nc.variables['uReconstructZonal'][0, :, :] = u_zonal_val
        nc.variables['uReconstructMeridional'][0, :, :] = u_merid_val
        nc.variables['pressure_p'][0, :, :] = pressure_val
        nc.variables['theta'][0, :, :] = theta_val
        nc.variables['qv'][0, :, :] = qv_val

# === Interpolate to 2D Lambert Grid ===
grid_surface_pressure = interp_irregular_to_regular(surface_pressure)
grid_uReconstructZonal = interp_3d_to_regular(uReconstructZonal)
grid_uReconstructMeridional = interp_3d_to_regular(uReconstructMeridional)
grid_pressure_p = interp_3d_to_regular(pressure_p)
grid_theta = interp_3d_to_regular(theta)
grid_qv = interp_3d_to_regular(qv)
# === Plot Example ===
def plot_field(grid, title, filename):
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(x_2d / 1000, y_2d / 1000, grid, shading='auto', cmap='viridis')
    plt.colorbar(label=title)
    plt.title(f"{title} at level {level_index}")
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

plot_field(grid_surface_pressure, "Surface Pressure", "surface_pressure_lambert.png")
plot_field(grid_uReconstructZonal[level_index], "Zonal Wind", "uReconstructZonal_lambert.png")
plot_field(grid_uReconstructMeridional[level_index], "Meridional Wind", "uReconstructMeridional_lambert.png")
plot_field(grid_pressure_p[level_index], "Pressure", "pressure_p_lambert.png")
plot_field(grid_theta[level_index], 'Theta', 'theta_lambert.png')
plot_field(grid_qv[level_index], 'QV', 'qv_lambert.png')

# === Write interpolated grid to NetCDF ===
with Dataset(output_nc_file, 'w') as nc_out:
    nc_out.createDimension('x', nx)
    nc_out.createDimension('y', ny)
    nc_out.createDimension('nVertLevels', nVertLevels)

    x_var = nc_out.createVariable('x', 'f4', ('x',))
    y_var = nc_out.createVariable('y', 'f4', ('y',))
    lev_var = nc_out.createVariable('nVertLevels', 'i4', ('nVertLevels',))
    x_var[:] = x_grid
    y_var[:] = y_grid
    lev_var[:] = np.arange(nVertLevels)

    x_var.units = 'meters'
    y_var.units = 'meters'
    x_var.standard_name = 'projection_x_coordinate'
    y_var.standard_name = 'projection_y_coordinate'

    # Add projection metadata
    proj = nc_out.createVariable('lambert_conformal_conic', 'i4')
    proj.grid_mapping_name = "lambert_conformal_conic"
    proj.standard_parallel = [lat0]
    proj.latitude_of_projection_origin = lat0
    proj.longitude_of_central_meridian = lon0
    proj.false_easting = 0.0
    proj.false_northing = 0.0
    proj.semi_major_axis = 6378137.0
    proj.inverse_flattening = 298.257223563

    def write_var(name, data, dims, units, long_name):
        var = nc_out.createVariable(name, 'f4', dims, zlib=True)
        var[:] = data
        var.units = units
        var.long_name = long_name
        var.grid_mapping = 'lambert_conformal_conic'

    write_var('surface_pressure', grid_surface_pressure, ('y', 'x'), 'Pa', 'Surface pressure')
    write_var('uReconstructZonal', grid_uReconstructZonal, ('nVertLevels', 'y', 'x'), 'm s-1', 'Zonal wind')
    write_var('uReconstructMeridional', grid_uReconstructMeridional, ('nVertLevels', 'y', 'x'), 'm s-1', 'Meridional wind')
    write_var('pressure_p', grid_pressure_p, ('nVertLevels', 'y', 'x'), 'Pa', 'Pressure')
    write_var('theta', grid_theta, ('nVertLevels', 'y', 'x'), 'K', 'Potential temperature')

print(f"✅ Lambert grid NetCDF written: {output_nc_file}")

# === Replace original MPAS variables with interpolated-back values ===
print(f"🔁 Creating {mpas_copy_a} and replacing variables with back-interpolated values...")


interp_surface = interp_regular_to_irregular(grid_surface_pressure)
interp_u = interp_3d_back(grid_uReconstructZonal)
interp_v = interp_3d_back(grid_uReconstructMeridional)
interp_p = interp_3d_back(grid_pressure_p)
interp_theta = interp_3d_back(grid_theta)
interp_qv = interp_3d_back(grid_qv)

# MPAS_COPY_A: Store smoothed data
update_mpas_file(
    mpas_copy_a,
    interp_surface+surface_pressure_mean,
    interp_u+uReconstructZonal_mean,
    interp_v+uReconstructMeridional_mean,
    interp_p+pressure_p_mean,
    interp_theta+theta_mean,
    interp_qv+qv_mean
)
print(f"✅ Variables updated in: {mpas_copy_a}")
# MPAS_COPY_B: Store residuals 
update_mpas_file(
    mpas_copy_b,
    surface_pressure-interp_surface+surface_pressure_mean,
    uReconstructZonal-interp_u+uReconstructZonal_mean,
    uReconstructMeridional-interp_v+uReconstructMeridional_mean,
    pressure_p-interp_p+pressure_p_mean,
    theta-interp_theta+theta_mean,
    qv-interp_qv+qv_mean
)
print(f"✅ Variables updated in: {mpas_copy_b}")

