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

mkdir -p togrid 
# Loop through all input files matching the pattern YYYYMMDD_HH.csv
#ioda.viirs_n20.2020093023.nc4
for input_file in ioda.viirs_npp.202009*.nc4; do
    # Extract the YYYYMMDDHH part from the filename
    datetime="${input_file:15:10}"  # Start at index 11, take 10 characters
    echo "$datetime"

    #mkdir -p $datetime
    ln -sf $input_file aod.nc 
    # Run the command with the input and output file names
    python3 iodaaod2dyngrid_yaml_v12.py >& togrid/log.${datetime}
    mv v12_interpolated_aod.nc  togrid/aod_$datetime.nc
    mv *.png togrid/aod_$datetime.png 
done

