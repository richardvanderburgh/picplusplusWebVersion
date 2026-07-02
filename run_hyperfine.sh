#!/usr/bin/env bash
set -euo pipefail

# Usage: ./run_hyperfine.sh [input_dir]
# Example: ./run_hyperfine.sh inputFiles/benchmarkFiles/particleDOE

input_dir="${1:-inputFiles/benchmarkFiles/timestepDOE}"
output_dir="results"
binary="${PIC_BINARY:-build/bin/PIC++Main}"

if [[ ! -x "$binary" && ! -f "$binary" ]]; then
  echo "Binary not found: $binary" >&2
  echo "Build first with ./scripts/build.sh or set PIC_BINARY." >&2
  exit 1
fi

mkdir -p "$output_dir"

for input_file in "$input_dir"/*.json; do
  base_name=$(basename "$input_file" .json)
  result_file="${output_dir}/${base_name}_result.json"

  hyperfine "$binary $input_file" --export-json "$result_file"
  echo "Completed benchmark for $input_file. Results saved to $result_file."
done

echo "All benchmarks completed."
