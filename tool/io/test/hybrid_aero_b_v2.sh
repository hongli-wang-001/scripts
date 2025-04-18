#!/bin/bash
set -ux 

# Define the input files and output file
file1="file1.nc"
file2="file2.nc"
output="output.nc"

# Create a temporary file to store the command script
script="ncap2_commands.nc"

# Get the list of variable names from file1 (assuming all variables are the same in both files)
variables=$(ncdump -h $file1 | grep -oP 'double\s+\K\w+' | tr '\n' ' ')

# Create the ncap2 command string
cmd="ncap2 -s "

for var in $variables; do
  cmd+="${var}=${var} + ${var}; "
done

# Remove the trailing semicolon and add input and output files
cmd=$(echo "$cmd" | sed 's/; $//')
cmd="$cmd $file1 $file2 $output"

# Execute the ncap2 command
eval $cmd

