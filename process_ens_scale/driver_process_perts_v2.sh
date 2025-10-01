#!/bin/bash

srcfile=/scratch3/NCEPDEV/fv3-cam/Ting.Lei/dr-jedi-bundle-current/src-jedi-bundle.sh
source $srcfile

# Set the ensemble directory and Python script
ENSEMBLE_DIR="data/ens"
PYTHON_SCRIPT="inp_lam_dct_v14.py"

mkdir -p data_ens/ens_a
mkdir -p data_ens/ens_b

# Loop over ensemble members 001 to 030
#for i in $(seq -w 1 1); do
for i in $(seq -f "%03g" 1 30); do
    echo "🔁 Processing ensemble member: mem${i}"

    # Construct the path to the forecast file
    FORECAST_FILE="${ENSEMBLE_DIR}/mem${i}/mpasout.2024-05-27_00.00.00.nc"

    # Check if forecast file exists
    if [ ! -f "$FORECAST_FILE" ]; then
        echo "⚠️ File not found: $FORECAST_FILE"
        continue
    fi

    # Remove existing symlink or file
    [ -e mpasout.nc ] && rm mpasout.nc

    # Create symbolic link
    ln -s "$FORECAST_FILE" mpasout.nc

    # Run the Python script
    echo "🚀 Running: python $PYTHON_SCRIPT"
    python "$PYTHON_SCRIPT"

    # Ensure target directories exist
    mkdir -p data_ens/ens_a/mem${i}
    mkdir -p data_ens/ens_b/mem${i}
    mkdir -p data_ens/ens_a/mem${i}/figure
    mkdir -p data_ens/ens_a/mem${i}/lambert
    # Move outputs
    mv mpasout_scale_a.nc data_ens/ens_a/mem${i}/mpasout.2024-05-27_00.00.00.nc
    mv mpasout_scale_b.nc data_ens/ens_b/mem${i}/mpasout.2024-05-27_00.00.00.nc
    mv interpolated_lambert.nc data_ens/ens_a/mem${i}/lam_mpasout.2024-05-27_00.00.00.nc
    mv *lambert.png data_ens/ens_a/mem${i}/lambert 
    mv *.png data_ens/ens_a/mem${i}/figure

    echo "✅ Finished mem${i}"
    echo "-----------------------------"
done

echo "🎉 All ensemble members processed."

