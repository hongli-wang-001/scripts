import h5py  # or netCDF4 if your file is in .nc format
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')

# Step 1: Load data from the NetCDF file
with h5py.File('file1.nc', 'r') as f1:
    tlon = f1['MetaData/longitude'][:]  # Longitude
    tlat = f1['MetaData/latitude'][:]   # Latitude
    pm = f1['ObsValue/particulatematter2p5Surface'][:]   # Observed PM2.5
    pm_sim = f1['hofx0/particulatematter2p5Surface'][:]  # Simulated PM2.5

# Step 2: Flatten the 2D or 3D data arrays to 1D
pm_flat = pm.flatten()      # Observed PM2.5 (flattened)
pm_sim_flat = pm_sim.flatten()  # Simulated PM2.5 (flattened)

# Step 3: Store the original indices for reshaping after transformation
# We will use argsort to get the original indices before transformation
original_indices = np.argsort(pm_sim_flat)

# Step 4: Sort the data to prepare for quantile mapping
sorted_pm = np.sort(pm_flat)         # Sorted observed PM2.5
sorted_pm_sim = np.sort(pm_sim_flat) # Sorted simulated PM2.5

# Step 5: Calculate the CDFs for both datasets (normalized ranks)
cdf_pm = np.linspace(0, 1, len(sorted_pm))  # CDF for observed data
cdf_pm_sim = np.linspace(0, 1, len(sorted_pm_sim))  # CDF for simulated data

# Step 6: Map the CDF of simulated data to the CDF of observed data (Quantile Mapping)
pm_sim_transformed = np.interp(cdf_pm_sim, cdf_pm, sorted_pm)

# Step 7: Reorder the transformed simulated data to match the original spatial grid
# We use the stored indices (original_indices) to reorder the transformed data
pm_sim_transformed_reordered = np.empty_like(pm_sim_flat)
pm_sim_transformed_reordered[original_indices] = pm_sim_transformed

# Step 8: Reshape the reordered transformed data back to the original shape of pm_sim
pm_sim_transformed_reordered = pm_sim_transformed_reordered.reshape(pm_sim.shape)

# Step 9: Visualize the original and transformed data distributions using histograms (PDF approximation)
plt.figure(figsize=(12, 6))

# Plot histogram of observed PM2.5 (PDF approximation)
plt.subplot(1, 2, 1)
plt.hist(pm_flat, bins=50, color='blue', alpha=0.7, label='Observed PM2.5', density=True)
plt.title('Histogram (PDF) of Observed PM2.5')
plt.xlabel('PM2.5 Concentration')
plt.ylabel('Probability Density')

# Plot histogram of transformed simulated PM2.5 (PDF approximation)
plt.subplot(1, 2, 2)
plt.hist(pm_sim_transformed_reordered.flatten(), bins=50, color='green', alpha=0.7, label='Transformed Simulated PM2.5', density=True)
plt.title('Histogram (PDF) of Transformed Simulated PM2.5')
plt.xlabel('PM2.5 Concentration')
plt.ylabel('Probability Density')

plt.tight_layout()
plt.savefig('q_mapping_pdf.png',dpi=300)
plt.show()

# Step 10: Visualize the spatial plot of transformed simulated data
myvalues=[2,3, 5, 7, 10, 15, 20, 30, 50, 70, 100]
plt.figure(figsize=(8, 6))
ax=plt.axes(projection=ccrs.PlateCarree())
plt.scatter(tlon,tlat,s=2,c=pm_sim_transformed_reordered,\
   linewidth=.50,cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
plt.axis([-135,-58,24,55])
cbar=plt.colorbar(fraction=0.07,orientation='horizontal')
cbar.set_label('Surface PM2.5 ($\mu$g/m$^3$)',fontsize=16)
plt.title('maped_pm',fontsize=18)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')

plt.savefig('q_mapping.png',dpi=300)

plt.show()

# Step 11: Optionally, save the transformed data to a new NetCDF file
# If you want to save the transformed data:
# with h5py.File('file_transformed.nc', 'w') as f_out:
#     f_out.create_dataset('Transformed_ParticulateMatter2p5Surface', data=pm_sim_transformed_reordered)

