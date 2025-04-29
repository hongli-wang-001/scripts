#!/bin/bash

#setenv PYTHONPATH "${PYTHONPATH}://scratch2/BMC/wrfruc/rli/JEDI/ioda-bundle/iodaconv/src:/scratch2/BMC/wrfruc/rli/JEDI/ioda-bundle/build/lib/python3.10"
#setenv sitefile False
module purge
module use /scratch1/NCEPDEV/jcsda/jedipara/spack-stack/modulefiles
module load miniconda/3.9.12
module load ecflow/5.11.4
module load mysql/8.0.31
module use /scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env-rocky8/install/modulefiles/Core
module load intel/2022.1.2
module load stack-intel/2021.5.0
module load impi/2022.1.2
module load stack-intel-oneapi-mpi/2021.5.1
module load stack-python/3.10.13
module load jedi-fv3-env
module load ewok-env

mkdir -p surfacepm2grid 
# Loop through all input files matching the pattern YYYYMMDD_HH.csv
for input_file in AIRNOW_202009????.nc; do
    # Extract the YYYYMMDDHH part from the filename
    datetime="${input_file:7:10}"  # Start at index 11, take 10 characters
    echo "$datetime"

    #mkdir -p $datetime
    ln -sf $input_file pm25.nc 
    # Run the command with the input and output file names
    python surfacebc2dyngrid_fill_v11.py >& surfacepm2grid/log.${datetime}
    mv v10_interpolated_pm25.nc  surfacepm2grid/bc_$datetime.nc
    mv *.png surfacepm2grid/bc_$datetime.png 
done

