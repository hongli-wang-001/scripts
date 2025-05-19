#!/bin/bash
set -x  # Echo commands for debugging

# Clean previous symlinks/files if they exist
[ -f fcst.nc ] && rm fcst.nc
[ -f pm25.nc ] && rm pm25.nc

# Configuration
export std=2020090200
export etd=2020090200
exp="3hcyc_pmf"

export fcstdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm_3hcyc/nco_dirs_${exp}/com/aqm/v7.0"
export obsdir="/scratch2/BMC/wrfruc/hwang/fire2_jedi/todyngrid_an/surfacepm2dyngrid"
export out="./vrf_${exp}"

mkdir -p "$out"

# Outer loop over dates
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
        # Format forecast hour as three-digit (e.g. 006)
        printf -v sth_padded "%03d" "$sth"

        fcstfile1="${fcstdir}/aqm.${year}${month}${dd}/${cyc}/aqm.t${cyc}z.dyn.f${sth_padded}.nc"

        # Advance time
        vrf_date=$(/home/Hongli.Wang/bin_theia/da_advance_time.exe "$std" "${sth}h")
        vrf_date_pre=$(/home/Hongli.Wang/bin_theia/da_advance_time.exe "$vrf_date" -1h)
        echo "vrf_date= $vrf_date"

        # Prepare input files
        cp "$obsdir/pm25_${vrf_date_pre}.nc" pm25.nc
        ln -sf "$fcstfile1" fcst.nc

        # Run verification script
        python vrf_gridobs_v4.py > "$outdir/log.${vrf_date}_f${sth_padded}.txt" 2>&1

        # Move output files
        mv pm25.nc "$outdir/pm25.${vrf_date}_f${sth_padded}.nc"
        [ -f pm25_comparison.pdf ] && mv pm25_comparison.pdf "$outdir/pm25.${vrf_date}_f${sth_padded}.pdf"

        # Increment forecast hour
        sth=$((sth + 1))
    done

    # Advance to next day/cycle (assuming 24-hour step)
    std=$(/home/Hongli.Wang/bin_theia/da_advance_time.exe "$std" 24h)
done
