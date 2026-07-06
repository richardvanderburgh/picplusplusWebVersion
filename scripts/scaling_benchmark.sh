#!/usr/bin/env bash
#
# OpenMP strong-scaling benchmark for PIC++.
#
# Runs the same fixed-size problem with an increasing number of OpenMP threads
# and records the wall-clock time of the core time-integration loop (the part
# that contains the parallel particle kernels). Results are written to
# results/scaling.json for plotting with scripts/plot_scaling.py.
#
# Usage:
#   scripts/scaling_benchmark.sh [input.json] [repeats]
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

INPUT="${1:-inputFiles/benchmarkFiles/scaling/large.json}"
REPEATS="${2:-3}"
OUTPUT="${SCALING_OUTPUT:-results/scaling.json}"

BIN="${PICPP_BIN:-build/bin/PIC++Main}"
if [[ ! -x "$BIN" ]]; then
  echo "Error: $BIN not found. Build first with ./scripts/build.sh" >&2
  exit 1
fi

if [[ ! -f "$INPUT" ]]; then
  echo "Error: input file '$INPUT' not found." >&2
  exit 1
fi

# Detect physical/logical core count for a sensible default thread ladder.
if [[ "$(uname -s)" == "Darwin" ]] && command -v sysctl &>/dev/null; then
  MAX_CORES="$(sysctl -n hw.ncpu)"
elif command -v nproc &>/dev/null; then
  MAX_CORES="$(nproc)"
else
  MAX_CORES=8
fi

# Build a thread ladder of powers of two up to MAX_CORES, always including
# MAX_CORES itself.
THREADS=()
t=1
while (( t < MAX_CORES )); do
  THREADS+=("$t")
  t=$(( t * 2 ))
done
THREADS+=("$MAX_CORES")

# Optional metadata for cross-platform comparison plots.
PLATFORM="${SCALING_PLATFORM:-$(uname -s | tr '[:upper:]' '[:lower:]')}"
COMPILER="${SCALING_COMPILER:-unknown}"
if [[ "$COMPILER" == "unknown" && -x "$BIN" ]]; then
  if strings "$BIN" 2>/dev/null | grep -q "GCC:"; then
    COMPILER="gcc"
  elif [[ "$(uname -s)" == "Darwin" ]]; then
    COMPILER="apple-clang"
  fi
fi
NPROC="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 0)"

echo "PIC++ OpenMP strong-scaling benchmark"
echo "  platform: $PLATFORM"
echo "  compiler: $COMPILER"
echo "  nproc:    $NPROC"
echo "  input:   $INPUT"
echo "  repeats: $REPEATS (reporting best/min)"
echo "  threads: ${THREADS[*]}"
echo ""

# Warn if this is a serial build (OpenMP not compiled in).
if "$BIN" "$INPUT" 2>/dev/null | grep -q "serial build"; then
  echo "WARNING: this is a SERIAL build (OpenMP not found at configure time)." >&2
  echo "         Thread counts will have no effect. Rebuild with libomp/GCC." >&2
  echo "" >&2
fi

mkdir -p results
run_times_us=()

for n in "${THREADS[@]}"; do
  best=""
  for ((r = 0; r < REPEATS; r++)); do
    out="$(OMP_NUM_THREADS="$n" "$BIN" "$INPUT" 2>/dev/null)"
    us="$(echo "$out" | sed -n 's/^Time loop took \([0-9][0-9]*\) micro secs$/\1/p')"
    if [[ -z "$us" ]]; then
      echo "Error: could not parse timing for $n threads." >&2
      exit 1
    fi
    if [[ -z "$best" || "$us" -lt "$best" ]]; then
      best="$us"
    fi
  done
  run_times_us+=("$best")
  printf "  %3d threads: %10d us\n" "$n" "$best"
done

# Emit results/scaling.json
{
  echo "{"
  printf '  "platform": "%s",\n' "$PLATFORM"
  printf '  "compiler": "%s",\n' "$COMPILER"
  printf '  "nproc": %s,\n' "$NPROC"
  printf '  "input": "%s",\n' "$INPUT"
  printf '  "repeats": %s,\n' "$REPEATS"
  printf '  "threads": ['
  for i in "${!THREADS[@]}"; do
    printf '%s' "${THREADS[$i]}"
    (( i < ${#THREADS[@]} - 1 )) && printf ', '
  done
  printf '],\n'
  printf '  "time_loop_us": ['
  for i in "${!run_times_us[@]}"; do
    printf '%s' "${run_times_us[$i]}"
    (( i < ${#run_times_us[@]} - 1 )) && printf ', '
  done
  printf ']\n'
  echo "}"
} > "$OUTPUT"

echo ""
echo "Wrote $OUTPUT"
echo "Plot with: .venv/bin/python scripts/plot_scaling.py"
