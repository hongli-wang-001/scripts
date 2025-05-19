#!/bin/bash
set -x  # Echo commands for debugging

# Clean previous symlinks/files if they exist
[ -f fcst.nc ] && rm fcst.nc
[ -f pm25.nc ] && rm pm25.nc

# Configuration
export std=2020090200  # Start time in YYYYMMDDHH
export etd=2020090200  # End time in YYYYMMDDHH
exp="3hcyc_pmf"

export fcstdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm_3hcyc/nco_dirs_${exp}/com/aqm/v7.0"
export obsdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/todyngrid_an/surfacepm2dyngrid"
export out="./vrf_${exp}"

mkdir -p "$out"

# Outer loop over forecast initialization times
while [ "$std" -le "$etd" ]; do
    export outdir="${out}/${std}"
    mkdir -p "$outdir"
    
    sth=1
    eth=24

    year=${std:0:4}
    month=${std:4:2}
    dd=${std:6:2}
    cyc=${std:8:2}

    while [ "$sth" -le "$eth" ]; do
        # Format forecast hour as 3 digits (e.g., 006)
        printf -v sth_padded "%03d" "$sth"

        fcstfile1="${fcstdir}/aqm.${year}${month}${dd}/${cyc}/aqm.t${cyc}z.dyn.f${sth_padded}.nc"

        # Advance verification time
        vrf_date=$(date -d "${std:0:8} ${std:8:2} +${sth} hours" +"%Y%m%d%H")
        vrf_date_pre=$(date -d "${vrf_date:0:8} ${vrf_date:8:2} -1 hour" +"%Y%m%d%H")
        echo "vrf_date = $vrf_date"

        # Prepare input files
        cp "$obsdir/pm25_${vrf_date_pre}.nc" pm25.nc
        ln -sf "$fcstfile1" fcst.nc

        # Run Python verification script and log output
        python vrf_gridobs_v4.py > "$outdir/log.${vrf_date}_f${sth_padded}.txt" 2>&1

        # Save outputs
        mv pm25.nc "$outdir/pm25.${vrf_date}_f${sth_padded}.nc"
        [ -f pm25_comparison.pdf ] && mv pm25_comparison.pdf "$outdir/pm25.${vrf_date}_f${sth_padded}.pdf"

        # Increment forecast hour
        sth=$((sth + 1))
    done

    # Advance initialization time by 24 hours
    std=$(date -d "${std:0:8} ${std:8:2} +24 hours" +"%Y%m%d%H")
done
