#!/bin/bash

setenv PYTHONPATH "${PYTHONPATH}://scratch2/BMC/wrfruc/rli/JEDI/ioda-bundle/iodaconv/src:/scratch2/BMC/wrfruc/rli/JEDI/ioda-bundle/build/lib/python3.10"
setenv sitefile False

mkdir -p surfacepm2grid 
# Loop through all input files matching the pattern YYYYMMDD_HH.csv
for input_file in pa-v4-*.nc; do
    # Extract the YYYYMMDDHH part from the filename
    datetime="${input_file:6:10}"  # Start at index 11, take 10 characters
    echo "$datetime"
    # Extract the base name (without extension)
    base_name="${input_file%.csv}"  # e.g., 20200914_12

    mkdir -p $datetime
    ln -sf $input_file pm25.nc 
    # Run the command with the input and output file names
    python3 surfacepm2grid_fill_v10.py >& surfacepm2grid/log.${datetime}
    mv v10_interpolated_pm25.nc  surfacepm2grid/pm25_$datetime.nc
    mv *.png surfacepm2grid/pm25_$datetime.png 
done

