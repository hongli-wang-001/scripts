#!/bin/bash

setenv PYTHONPATH "${PYTHONPATH}://scratch2/BMC/wrfruc/rli/JEDI/ioda-bundle/iodaconv/src:/scratch2/BMC/wrfruc/rli/JEDI/ioda-bundle/build/lib/python3.10"
setenv sitefile False

mkdir -p ioda_out_pa_v4
# Loop through all input files matching the pattern YYYYMMDD_HH.csv
for input_file in *.csv; do
    # Extract the base name (without extension)
    base_name="${input_file%.csv}"  # e.g., 20200914_12

    # Construct the output file name
    output_file="pa-v4-${base_name//_/}.nc"  # e.g., pa-v3-2020091412.nc

    # Run the command with the input and output file names
    #python3 purpleair2ioda-nc-v3.py -i "$input_file" -o ioda_out_pa_v3/"$output_file" 
    python3 purpleair2ioda-nc-superob_v4.py -i "$input_file" -o ioda_out_pa_v4/"$output_file"
done

