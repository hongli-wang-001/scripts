#!/bin/csh

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
module load soca-env
#module list

foreach file_3dvar (3dvar_nicas_pm25_202009????.nc)
    # Create a symbolic link to the current file
    ln -sf ./$file_3dvar ./3dvar_nicas_pm25.nc

python3 plot_scatter_omb_oma_v9_1file.py >! log.$file_3dvar
python3 plot_frq_omb_oma_v4a.py
python3 obserr-pm25.py $file_3dvar OBS_ERROR OBSERR-PM2p5
python3 simulated-amb-pm25.py $file_3dvar A-B AMB-PM2p5
python3 obs-pm25_v2.py     $file_3dvar AirNOW_PM2.5 OBS_PM2p5
python3 simulated-bkg-pm25.py  $file_3dvar bkg_PM2.5 BKG_PM2p5
python3 simulated-ana-pm25.py  $file_3dvar ana_PM2.5 ANA_PM2p5
python3 omb_3dvar_pm25_v2.py $file_3dvar O-B OMB_PM2p5_VADER
python3 oma_3dvar_pm25_v2.py $file_3dvar O-A OMA_PM2p5_VADER

    # Create a directory for the output PNG files
    mkdir -p ./png_$file_3dvar

    # Move the PNG files into the created directory
    mv ./*.png png_$file_3dvar
    mv ./frequency_info.nc png_$file_3dvar
end

