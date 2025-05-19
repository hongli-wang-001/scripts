#!/bin/bash

# Template Python script
template="compute_2d_stats_v3_template.py"

# Check that template exists
if [[ ! -f $template ]]; then
    echo "Template script $template not found!"
    exit 1
fi

# Loop over experiments and hours
#for exp in vrf_3hcyc_pm_h100v5; do
#for exp in vrf_6hcyc_ctr vrf_6hcyc_pm_h100v5 vrf_6hcyc_pmaod_h100v5  ; do
#for exp in vrf_3hcyc_pa_h100v12 ; do
for exp in vrf_6hcyc_ctr  ; do
  for hh in 00 06 12 18 all; do  # use 'all' instead of '??' for clarity
    echo "Processing experiment $exp, cycle hour $hh..."

    # Create customized Python script
    script_temp="compute_2d_stats_v3_temp.py"
    script="compute_2d_stats_v3_${exp}_${hh}.py"

    if [[ "$hh" == "all" ]]; then
      sed "s/HH/??/g" "$template" > "$script_temp"
      sed "s/EXP/${exp}/g" "$script_temp" > "$script"
    else
      sed "s/HH/${hh}/g" "$template" > "$script_temp"
      sed "s/EXP/${exp}/g" "$script_temp" > "$script"
    fi 

    # Set output directory
    if [[ "$hh" == "all" ]]; then
      outdir="${exp}/iall"
    else
      outdir="${exp}/i${hh}"
    fi

    mkdir -p "$outdir"

    # Run the script
    #python "$script"
    python "$script" > "${outdir}/run.compute_2d_stats_v3_${exp}_${hh}.log" 2>&1

    if [[ $? -ne 0 ]]; then
      echo "Python script failed for $hh in $exp"
      #exit 2
    fi

    mv aod_2d*.nc "$outdir"

    # Clean up
    rm -f "$script_temp" #"$script"

    echo "Finished cycle $hh."
  done
done

