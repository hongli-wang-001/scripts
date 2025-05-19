#!/bin/bash
set -x  # Echo commands for debugging

# Clean previous symlinks/files if they exist
[ -f fcst.nc ] && rm fcst.nc
[ -f aod.nc ] && rm aod.nc

# Configuration
export std=2020090200  # Start time in YYYYMMDDHH
export etd=2020092118  # End time in YYYYMMDDHH
#exp="3hcyc_pmf"
#exp_name="3hcyc_pm_h100v12"
#exp="3hcyc_ctr"
#exp_name="3hcyc_ctr"
#exp="3hcyc_pmfvt5_airnow"
#exp_name="3hcyc_pm_h100v5"
exp="3hcyc_pmaod"
exp_name="3hcyc_pmaod_h100v5"
#exp="6hcyc_pmfvt5"
#exp_name="6hcyc_pm_h100v5"
#exp="6hcyc_ctr"
#exp_name="6hcyc_ctr"
#exp="6hcyc_pmaodvt5"
#exp_name="6hcyc_pmaod_h100v5"
#exp="3hcyc_pmfvt5_purpleair"
#exp_name="3hcyc_pa_h100vt5"
#exp="3hcyc_pmg"
#exp_name="3hcyc_pa_h100v12"
#exp="3hcyc_pmvPBL"
#exp_name="3hcyc_pm_h100vPBL_b1"

export fcstdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm_3hcyc/nco_dirs_${exp}/com/aqm/v7.0"
export obsdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/todyngrid_an/aeronet2dyngrid/aeroaod2grid"
export out="./vrf_${exp_name}"

mkdir -p "$out"

# Outer loop over forecast initialization times
while [ "$std" -le "$etd" ]; do
    export outdir="${out}/${std}"
    mkdir -p "$outdir"
    
    sth=0
    eth=23

    year=${std:0:4}
    month=${std:4:2}
    dd=${std:6:2}
    cyc=${std:8:2}

    while [ "$sth" -le "$eth" ]; do
        # Format forecast hour as 3 digits (e.g., 006)
        printf -v sth_padded "%03d" "$sth"

        fcstfile1="${fcstdir}/aqm.${year}${month}${dd}/${cyc}/aqm.t${cyc}z.phy.f${sth_padded}.nc"

        # Advance verification time
        vrf_date=$(date -d "${std:0:8} ${std:8:2} +${sth} hours" +"%Y%m%d%H")
        vrf_date_pre=$(date -d "${vrf_date:0:8} ${vrf_date:8:2} -1 hour" +"%Y%m%d%H")
        echo "vrf_date = $vrf_date"

        # Prepare input files
        #cp "$obsdir/aod_${vrf_date_pre}.nc" aod.nc
        cp "$obsdir/aod_${vrf_date}.nc" aod.nc
        ln -sf "$fcstfile1" fcst.nc

        # Run Python verification script and log output
        python vrf_gridaod_v6.py > "$outdir/log.${vrf_date}_f${sth_padded}.txt" 2>&1

        # Save outputs
        mv aod.nc "$outdir/aod.${vrf_date}_f${sth_padded}.nc"
        [ -f aod_comparison.pdf ] && mv aod_comparison.pdf "$outdir/aod.${vrf_date}_f${sth_padded}.pdf"

        # Increment forecast hour
        sth=$((sth + 1))
    done

    # Advance initialization time by 6 hours
    std=$(date -d "${std:0:8} ${std:8:2} +6 hours" +"%Y%m%d%H")
done
