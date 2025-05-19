#!/bin/bash
set -x  # Echo commands for debugging

# Enable `module` command in sh
source /apps/lmod/lmod/init/sh  # adjust path if needed

# Now modules will work
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

# Configuration
export std=2020090200  # Start time
export etd=2020092118  # End time
#export obsdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/todyngrid_an/iodaaod2grid/ioda.viirs.n20togrid"
export obsdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/todyngrid_an/surfacebc2dyngrid"
# /scratch2/BMC/wrfruc/hwang/fire2_jedi/todyngrid_an/surfacebc2dyngrid/bc_2020091418.nc

# Declare experiments and names as parallel arrays
experiments=(
  "3hcyc_pmf"
  "3hcyc_ctr"
  "3hcyc_pmfvt5_airnow"
  "3hcyc_pmaod"
  "6hcyc_pmfvt5"
  "6hcyc_ctr"
  "6hcyc_pmaodvt5"
  "3hcyc_pmfvt5_purpleair"
  "3hcyc_pmg"
  "3hcyc_pmvPBL"
)

exp_names=(
  "3hcyc_pm_h100v12"
  "3hcyc_ctr"
  "3hcyc_pm_h100v5"
  "3hcyc_pmaod_h100v5"
  "6hcyc_pm_h100v5"
  "6hcyc_ctr"
  "6hcyc_pmaod_h100v5"
  "3hcyc_pa_h100vt5"
  "3hcyc_pa_h100v12"
  "3hcyc_pm_h100vPBL_b1"
)

rm bc.nc
# Loop over each experiment
for i in "${!experiments[@]}"; do
  exp="${experiments[i]}"
  exp_name="${exp_names[i]}"
  fcstdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm_3hcyc/nco_dirs_${exp}/com/aqm/v7.0"
  out="./vrf_${exp_name}"
  mkdir -p "$out"

  curr_std="$std"
  while [ "$curr_std" -le "$etd" ]; do
    outdir="${out}/${curr_std}"
    mkdir -p "$outdir"

    sth=0
    eth=24

    year=${curr_std:0:4}
    month=${curr_std:4:2}
    dd=${curr_std:6:2}
    cyc=${curr_std:8:2}

    while [ "$sth" -le "$eth" ]; do
      printf -v sth_padded "%03d" "$sth"

      fcstfile1="${fcstdir}/aqm.${year}${month}${dd}/${cyc}/aqm.t${cyc}z.dyn.f${sth_padded}.nc"
      vrf_date=$(date -d "${curr_std:0:8} ${curr_std:8:2} +${sth} hours" +"%Y%m%d%H")
      vrf_date_pre=$(date -d "${vrf_date:0:8} ${vrf_date:8:2} -1 hour" +"%Y%m%d%H")

      echo "vrf_date = $vrf_date"

      # Check and copy observation file
      if [ -f "$obsdir/bc_${vrf_date}.nc" ]; then
        cp "$obsdir/bc_${vrf_date}.nc" bc.nc
        ln -sf "$fcstfile1" fcst.nc

        python vrf_gridbc_v6.py > "$outdir/log.${vrf_date}_f${sth_padded}.txt" 2>&1

        mv bc.nc "$outdir/bc.${vrf_date}_f${sth_padded}.nc"
        [ -f bc_comparison.pdf ] && mv bc_comparison.pdf "$outdir/bc.${vrf_date}_f${sth_padded}.pdf"
      else
        echo "Warning: File $obsdir/bc_${vrf_date}.nc does not exist."
      fi

      sth=$((sth + 3))
    done

    curr_std=$(date -d "${curr_std:0:8} ${curr_std:8:2} +6 hours" +"%Y%m%d%H")
  done
done

