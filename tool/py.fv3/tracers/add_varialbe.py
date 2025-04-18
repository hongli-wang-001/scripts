import netCDF4 as nc
import numpy as np

# Define the input and output file names
file1 = 'file1.nc'
file2 = 'file2.nc'
output_file = 'output.nc'

# Open the input NetCDF files
with nc.Dataset(file1, 'r') as ds1, nc.Dataset(file2, 'r') as ds2:
    # Check that both files have the same variables
    variables = list(ds1.variables.keys())
    
    # Ensure that both files contain the same variables
    if variables != list(ds2.variables.keys()):
        raise ValueError("The two files do not contain the same variables.")
    
    # Create a new NetCDF file for output
    with nc.Dataset(output_file, 'w', format='NETCDF4') as ds_out:
        # Copy dimensions from the first file
        for name, dimension in ds1.dimensions.items():
            ds_out.createDimension(name, len(dimension) if not dimension.isunlimited() else None)
        
        # Copy variables and add them
        for var_name in variables:
            var1 = ds1.variables[var_name][:]
            var2 = ds2.variables[var_name][:]
            
            # Check if the dimensions match
            if var1.shape != var2.shape:
                raise ValueError(f"Dimensions of variable '{var_name}' do not match between files.")
            
            # Create the new variable in the output file
            var_out = ds_out.createVariable(f'{var_name}_sum', var1.dtype, ds1.variables[var_name].dimensions)
            
            # Add the variables from the two files
            var_out[:] = var1 + var2

            # Copy attributes
            var_out.setncattr('long_name', f'Sum of {var_name} from file1 and file2')
            for attr_name in ds1.variables[var_name].ncattrs():
                var_out.setncattr(attr_name, ds1.variables[var_name].getncattr(attr_name))

print(f"Finished adding variables and saved to {output_file}")

