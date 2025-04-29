#!/bin/bash

# Define start and end cycle timestamps
start_cycle="202009110600"
end_cycle="202009211800"

# Path to your workflow XML and DB
workflow_xml="FV3LAM_wflow.xml"
workflow_db="FV3LAM_wflow.xml.db"

# Convert to seconds since epoch for arithmetic
start_sec=$(date -d "${start_cycle:0:4}-${start_cycle:4:2}-${start_cycle:6:2} ${start_cycle:8:2}:00:00" +%s)
end_sec=$(date -d "${end_cycle:0:4}-${end_cycle:4:2}-${end_cycle:6:2} ${end_cycle:8:2}:00:00" +%s)

# Loop in 3-hour steps (10800 seconds)
step_sec=10800

while [ "$start_sec" -le "$end_sec" ]; do
  cycle=$(date -u -d "@$start_sec" +%Y%m%d%H%M)

  echo "Checking cycle: $cycle"
  rocotocheck -v10 -w "$workflow_xml" -d "$workflow_db" -c "$cycle" -t gen_yaml
  echo "------------------------------------------------------"

  start_sec=$((start_sec + step_sec))
done
