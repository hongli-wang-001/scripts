#!/bin/bash

# Define the input files and output file
file1="file1.nc"
file2="file2.nc"
output="output.nc"

# Start with an empty command
cmd="ncap2 -s "

# Loop over the variables
for var in $(ncdump -h $file1 | grep "double" | awk '{print $2}' | sed 's/;//'); do
  cmd+="\"${var}_sum=${var} + ${var}\"\n"
done

# Remove the trailing newline
cmd=$(echo "$cmd" | sed 's/;$//')

# Add the filenames to the command
cmd+=" $file1 $file2 $output"

# Execute the command
eval $cmd
