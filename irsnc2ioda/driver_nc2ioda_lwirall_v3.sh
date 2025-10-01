#!/bin/bash

set -euo pipefail

# Usage/help function
usage(){
  echo "Usage: $0 [--skip-step1] [--skip-step2]"
  echo
  echo "Options:"
  echo "  --skip-step1   Skip Step 1 (processing input files / generating individual IODA files)"
  echo "  --skip-step2   Skip Step 2 (making Location dimension unlimited on IODA files)"
  exit 1
}

# Default: run both steps
DO_STEP1=true
DO_STEP2=false

# Parse flags
while (( "$#" )); do
  case "$1" in
    --skip-step1)
      DO_STEP1=false
      shift
      ;;
    --skip-step2)
      DO_STEP2=false
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
done

INPUT_DIR="input_5"
PATTERN="*1B-PC*BODY*.nc"
PY_SCRIPT="nc2ioda_derived_rad_v29.py"
LINK_NAME="derived_rad.nc"
OUTPUT_BASE="derived_rad_v29.ioda.nc"
IODA_PATTERN="ioda_DEV_*.ioda.nc"
FINAL_OUTPUT="combined_ioda.nc"

# Step 1: produce individual IODA files
if [ "$DO_STEP1" = true ]; then
  echo "=== Step 1: Generating individual IODA files ==="
  for input_file in "$INPUT_DIR"/$PATTERN; do
    [ -f "$input_file" ] || continue
    echo "Processing: $input_file"
    rm -f "$LINK_NAME" "$OUTPUT_BASE"
    ln -s "$(realpath "$input_file")" "$LINK_NAME"
    python3 "$PY_SCRIPT"
    if [ ! -f "$OUTPUT_BASE" ]; then
      echo "ERROR: Output $OUTPUT_BASE not found for $input_file"
      continue
    fi
    # Extract times from filename
    fname=$(basename "$input_file")
    if [[ $fname =~ _([0-9]{14})_([0-9]{14})_ ]]; then
      tstart="${BASH_REMATCH[1]}"
      tend="${BASH_REMATCH[2]}"
    else
      tstart="UNKNOWNSTART"
      tend="UNKNOWNEND"
    fi
    newname="ioda_DEV_${tstart}_${tend}.ioda.nc"
    mv "$OUTPUT_BASE" "$newname"
    echo "Renamed output to $newname"
  done
else
  echo "=== Skipping Step 1 ==="
fi

# Step 2: ensure Location dimension is unlimited
if [ "$DO_STEP2" = true ]; then
  echo "=== Step 2: Making Location unlimited where needed ==="
  for f in $IODA_PATTERN; do
    [ -f "$f" ] || continue
    # Check if Location dimension is already unlimited
    if ncdump -h "$f" | grep -q "Location = UNLIMITED" ; then
      echo "File $f: Location already unlimited. Skipping."
    else
      echo "File $f: Making Location unlimited."
      ncks --mk_rec_dmn Location "$f" tmp_"$f"
      mv tmp_"$f" "$f"
    fi
  done
else
  echo "=== Skipping Step 2 ==="
fi

# Step 3: concatenate into one
echo "=== Step 3: Concatenating all IODA files ==="
ncrcat $IODA_PATTERN $FINAL_OUTPUT
echo "Concatenation done: $FINAL_OUTPUT"


