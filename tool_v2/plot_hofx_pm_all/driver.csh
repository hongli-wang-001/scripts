#!/bin/csh
module use /scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm/aqm_jedi/modulefiles
module use /contrib/miniconda3/modulefiles 
module load wflow_hera
conda activate regional_workflow

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
module unload fms
module load fms/2023.04


foreach dir ( 3dvar_dif*_pm25_2020092*.nc)
ln -sf $dir file1.nc 
python3 pdf_log1ptrans_kdefit_v4b.py
mkdir -p figure_hofx_$dir
mv *.png figure_hofx_$dir 
driver_plot_3dvar.csh $dir
end

