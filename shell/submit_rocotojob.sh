#!/bin/bash

start_cycle="202009110300"
end_cycle="202009211800"

workflow_xml="FV3LAM_wflow.xml"
workflow_db="FV3LAM_wflow.xml.db"

# Convert to seconds for loop
start_sec=$(date -d "${start_cycle:0:4}-${start_cycle:4:2}-${start_cycle:6:2} ${start_cycle:8:2}:00:00" +%s)
end_sec=$(date -d "${end_cycle:0:4}-${end_cycle:4:2}-${end_cycle:6:2} ${end_cycle:8:2}:00:00" +%s)

step_sec=10800  # 3 hours

while [ "$start_sec" -le "$end_sec" ]; do
  cycle=$(date -u -d "@$start_sec" +%Y%m%d%H%M)

  # Extract date/time components
  yyyy=${cycle:0:4}
  mm=${cycle:4:2}
  dd=${cycle:6:2}
  hh=${cycle:8:2}
  yyyymmdd=${yyyy}${mm}${dd}
  hh_cycle=$(printf "%02d" $((10#$hh)))
  model_time="${yyyy}${mm}${dd}.${hh}0000"

  prev_sec=$((start_sec - step_sec))
  prev_cycle=$(date -u -d "@$prev_sec" +%Y%m%d%H%M)

  # Extract date/time components
  yyyy2=${prev_cycle:0:4}
  mm2=${prev_cycle:4:2}
  dd2=${prev_cycle:6:2}
  hh2=${prev_cycle:8:2}
  prev_yyyymmdd=${yyyy2}${mm2}${dd2}
  

  # Expected model file
  model_file="/scratch2/BMC/wrfruc/hwang/fire2_jedi/aqm_3hcyc/expt_dirs/../nco_dirs_3hcyc_pmfvt5/com/aqm/v7.0/aqm.${prev_yyyymmdd}/${hh2}/RESTART/${model_time}.fv_tracer.res.tile1.nc"

  echo "Checking cycle: $cycle"
  
  # Run rocotocheck and capture output
  check_output=$(rocotocheck -v10 -w "$workflow_xml" -d "$workflow_db" -c "$cycle" -t gen_yaml 2>&1)

  # Check conditions
  if echo "$check_output" | grep -q "The cycle is not active"; then
    #if echo "$check_output" | grep -q "fv_tracer.res.tile1.nc is available"; then
    if [ -f "$model_file" ]; then
      echo "🟡 Cycle $cycle is inactive, but model file exists."
      echo "🔁 Rewinding cycle $cycle..."
      /apps/rocoto/1.3.7/bin/rocotoboot -w "$workflow_xml" -d "$workflow_db" -c "$cycle"  -t gen_yaml
      echo "✅ Exiting loop after first rewind."
      break
    else
      echo "🔴 Model file not found: $model_file"
    fi
  else
    echo "✅ Task is active or already submitted."
  fi

  echo "------------------------------------------------------"

  # Step forward
  start_sec=$((start_sec + step_sec))
done
/apps/rocoto/1.3.7/bin/rocotorun  -v10  -w FV3LAM_wflow.xml -d FV3LAM_wflow.xml.db
echo "✅ rocotorun submitted."
