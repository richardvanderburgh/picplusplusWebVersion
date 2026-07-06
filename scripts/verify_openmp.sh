#!/usr/bin/env bash
#
# Verify that PIC++Main was built with OpenMP enabled.
# Exits 0 when the binary reports an active OpenMP thread pool; exits 1 for a
# serial build or if the binary is missing.
#
# Usage:
#   scripts/verify_openmp.sh [binary] [input.json]
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

BIN="${1:-build/bin/PIC++Main}"
INPUT="${2:-inputFiles/benchmarkFiles/scaling/ci_smoke.json}"

if [[ ! -x "$BIN" ]]; then
  for candidate in \
    build/bin/PIC++Main \
    build/bin/PIC++Main.exe \
    build/bin/Release/PIC++Main.exe; do
    if [[ -x "$candidate" ]]; then
      BIN="$candidate"
      break
    fi
  done
fi

if [[ ! -x "$BIN" ]]; then
  echo "ERROR: executable not found: $BIN" >&2
  exit 1
fi

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: input file not found: $INPUT" >&2
  exit 1
fi

# Only parse the startup banner; avoid running a long simulation when possible.
out="$(OMP_NUM_THREADS=2 "$BIN" "$INPUT" 2>&1 | head -n 20)"

if echo "$out" | grep -q "serial build"; then
  echo "FAIL: $BIN is a serial build (OpenMP not enabled at compile time)." >&2
  if [[ "$(uname -s)" == "Darwin" ]]; then
    echo "      On macOS install libomp, then rebuild: brew install libomp && ./scripts/build.sh" >&2
  fi
  exit 1
fi

threads="$(echo "$out" | sed -n 's/^OpenMP threads: \([0-9][0-9]*\)$/\1/p' | head -n 1)"
if [[ -z "$threads" || "$threads" -lt 2 ]]; then
  echo "FAIL: expected OpenMP threads >= 2 with OMP_NUM_THREADS=2, got: '${threads:-<none>}'" >&2
  echo "$out" >&2
  exit 1
fi

echo "OK: OpenMP enabled ($threads threads reported with OMP_NUM_THREADS=2)"
exit 0
