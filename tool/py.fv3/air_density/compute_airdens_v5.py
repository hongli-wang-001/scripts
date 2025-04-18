import numpy as np
from netCDF4 import Dataset

def compute_air_density(temp, delp, rdgas=287.58):
    """
    Compute the air density from temperature and pressure difference.
    
    Parameters:
    - temp: 3D numpy array of temperature (lev, lat, lon)
    - delp: 3D numpy array of pressure difference (lev, lat, lon)
    - rdgas: Specific gas constant for dry air
    
    Returns:
    - airdens: 3D numpy array of air density (lev, lat, lon)
    """
    # Compute cumulative pressure levels
    pfull = np.cumsum(delp, axis=0)
    
    # Calculate air density
    airdens = pfull / (rdgas * temp)
    
    return airdens

def main():
    # Open the NetCDF files
    A = Dataset("./airdens.nc", "r+")  # File to write results
    B = Dataset("./gfsic.nc", "r")     # File containing temperature, pressure
    C = Dataset("./airdens_fv.nc", "r") # File containing reference air density

    # Read variables
    temp = B.variables['t'][:]
    delp = B.variables['delp'][:]
    ps = B.variables['ps'][:]
    airdens_r = C.variables['dry_air_density'][:]

    # Compute air density
    airdens = compute_air_density(temp, delp)

    # Convert airdens_r to float to match the type of airdens
    airdens_r = airdens_r.astype(np.float32)

    # Initialize diff
    diff = np.full_like(airdens, np.nan)

    # Calculate the difference
    diff[1:65, :, :] = airdens[1:65, :, :] - airdens_r[0:64, :, :]

    # Check if the variable already exists and delete it if so
    if 'diff_dry_air_density' in A.variables:
        print("Deleting existing variable: diff_dry_air_density")
        del A.variables['diff_dry_air_density']
    
    # Create the new variable
    diff_var = A.createVariable('diff_dry_air_density', 'f4', ('lev', 'lat', 'lon'))

    # Assign data to the new variable
    diff_var[:] = diff

    # Add attributes to the new variable (if needed)
    diff_var.long_name = 'Difference of dry air density'
    diff_var.units = 'kg m-3'

    # Close the files
    A.close()
    B.close()
    C.close()

if __name__ == "__main__":
    main()

