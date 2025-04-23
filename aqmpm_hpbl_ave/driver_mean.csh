#!/bin/csh -x
#
#
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
mkdir -p i00
mkdir -p i06
mkdir -p i12
mkdir -p i18
foreach datel (001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020 021 022 023 024) 
#sed s/f001/f$datel/g mean_std_of_diff_excf000.i18z.con.vs.pmf.v3.py > mean_std_of_diff_excf000.i18z.con.vs.pmf.v3_f$datel.py
#python3 mean_std_of_diff_excf000.i18z.con.vs.pmf.v3_f$datel.py
#mv output_statistics.nc i18/statistics.f$datel.nc

sed s/f001/f$datel/g mean_std_of_diff_excf000.i12z.con.vs.pmf.v3.py > mean_std_of_diff_excf000.i12z.con.vs.pmf.v3_f$datel.py
python3 mean_std_of_diff_excf000.i12z.con.vs.pmf.v3_f$datel.py
mv output_statistics.nc i12/statistics.f$datel.nc

sed s/f001/f$datel/g mean_std_of_diff_excf000.i06z.con.vs.pmf.v3.py > mean_std_of_diff_excf000.i06z.con.vs.pmf.v3_f$datel.py
python3 mean_std_of_diff_excf000.i06z.con.vs.pmf.v3_f$datel.py
mv output_statistics.nc i06/statistics.f$datel.nc

sed s/f001/f$datel/g mean_std_of_diff_excf000.i00z.con.vs.pmf.v3.py > mean_std_of_diff_excf000.i00z.con.vs.pmf.v3_f$datel.py
python3 mean_std_of_diff_excf000.i00z.con.vs.pmf.v3_f$datel.py
mv output_statistics.nc i00/statistics.f$datel.nc

end

