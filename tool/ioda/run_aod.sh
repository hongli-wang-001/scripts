#!/bin/bash -l
#
#
# -- Specify a maximum wallclock of 4 hours
#SBATCH --account=zrtrr
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=24000
#SBATCH -o log.aod

#set -ux

#source /scratch1/BMC/wrfruc/gge/skylab/6.0/modulefile.skylab6
#export PYTHONPATH=$PYTHONPATH:/scratch1/BMC/wrfruc/bjallen/skylab/build/lib/python3.10/
source /home/Ruifang.Li/modulefile.skylab6.rocky8

export SLURM_ACCOUNT=zrtrr
export SALLOC_ACCOUNT=$SLURM_ACCOUNT
export SBATCH_ACCOUNT=$SLURM_ACCOUNT


aod_j01=/scratch1/BMC/wrfruc/hwang/wf1/data/case2020sep/aod_viis/j01/001
aod_npp=/scratch1/BMC/wrfruc/hwang/wf1/data/case2020sep/aod_viis/npp/001

yy=2020
mm=09

for dd in  $(seq -w 13 13); do
for cyc in $(seq -w 18 18); do	

  for path in `ls ${aod_npp}/JRR-AOD_*s${yy}${mm}${dd}${cyc}*nc`; do
    file=$(basename "$path")
    #srun python3 -u ./viirs_aod2ioda.py  -i ${aod_j01}/JRR-AOD_v2r3_j01_s202009121815*.nc -o ./Data/output/JRR-AOD_v2r3_j01_s202009121815.nc -m nesdis -k maskout 
    srun python3 -u ./viirs_aod2ioda_org.py  -i ${aod_npp}/$file -o ./Data/output/$file -m nesdis -k maskout 
  done
done
done

# run from login node, need to setup env 
# python3  ./viirs_aod2ioda.py  --adp_mask  -i /scratch1/BMC/wrfruc/hwang/wf1/data/case2020sep/aod_viis/j01/001/JRR-AOD_v2r3_j01_s202009302259508*.nc -o ./Data/output/JRR-AOD_v2r3_j01_s20200930.nc -m nesdis -k maskout

# python3  ./viirs_aod2ioda.py  --adp_mask  -i /scratch1/BMC/wrfruc/bjallen/skylab/viirs/aod/conus/JRR-AOD_v3r2_npp_s2024010711*.nc -o ./Data/output/JRR-AOD_v2r3_npp_test.nc -m nesdis -k maskout
