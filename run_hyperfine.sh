#!/bin/bash

# Define the directory containing your input files
input_dir="inputFiles/benchmarkFiles/timestepDOE"

# Define the output directory for your results files
output_dir="results"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Get a list of input files
input_files=$(ls $input_dir/*.json)

# Iterate over each input file
for input_file in $input_files; do
    # Extract the basename of the file for naming the result file
    base_name=$(basename "$input_file" .json)
    
    # Define the output file name
    result_file="${output_dir}/${base_name}_result.json"
    
    # Run Hyperfine and export the results to a JSON file
    hyperfine "build/bin/PIC++Main $input_file" --export-json "$result_file"
    
    # Optional: Echo the completion of this test
    echo "Completed benchmark for $input_file. Results saved to $result_file."
done

echo "All benchmarks completed."
