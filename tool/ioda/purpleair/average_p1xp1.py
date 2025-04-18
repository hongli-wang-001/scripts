import numpy as np
import pandas as pd

# Assuming data is a structured array, convert it to a DataFrame for easier manipulation
df = pd.DataFrame(data)

# Define the grid size for averaging
lat_bins = np.arange(df['latitude'].min(), df['latitude'].max() + 0.1, 0.1)
lon_bins = np.arange(df['longitude'].min(), df['longitude'].max() + 0.1, 0.1)

# Create a new column for latitude and longitude bin labels
df['lat_bin'] = pd.cut(df['latitude'], bins=lat_bins)
df['lon_bin'] = pd.cut(df['longitude'], bins=lon_bins)

# Group by the bin labels and calculate the average for var1 and var2
averaged_data = df.groupby(['lat_bin', 'lon_bin'])[[var1, var2]].mean().reset_index()

# Optionally, rename the averaged columns if needed
averaged_data = averaged_data.rename(columns={var1: f'avg_{var1}', var2: f'avg_{var2}'})

# Convert the averaged data back to a structured array if needed
averaged_array = np.zeros(averaged_data.shape[0], dtype=data.dtype)
for col in averaged_data.columns:
    if col in data.dtype.names:
        averaged_array[col] = averaged_data[col].values

# Now you can use averaged_array for further processing or output

