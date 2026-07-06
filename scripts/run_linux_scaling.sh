#!/usr/bin/env bash
#
# Run a Linux GCC strong-scaling benchmark.
#
# Preferred path (no local Docker): trigger the GitHub Actions workflow and
# download the artifact:
#   gh workflow run scaling_linux.yml
#   gh run list --workflow=scaling_linux.yml --limit 1
#   gh run download <run-id> -n scaling-linux-json -D results/
#
# Local Docker fallback (when docker is installed):
#   ./scripts/run_linux_scaling.sh --docker
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [[ "${1:-}" == "--docker" ]]; then
  if ! command -v docker &>/dev/null; then
    echo "docker not found. Use: gh workflow run scaling_linux.yml" >&2
    exit 1
  fi
  docker run --rm --platform linux/amd64 \
    -v "$ROOT:/workspace" -w /workspace \
    ubuntu:22.04 bash -lc '
      set -euo pipefail
      apt-get update -qq
      DEBIAN_FRONTEND=noninteractive apt-get install -y -qq build-essential cmake python3 python3-venv git
      python3 -m venv .venv-linux
      .venv-linux/bin/pip install -q "conan>=2.0" cmake
      export PATH="$PWD/.venv-linux/bin:$PATH"
      PROFILE=buildUtils/linux_gcc_release ./scripts/build.sh
      SCALING_OUTPUT=results/scaling_linux.json \
      SCALING_PLATFORM=linux \
      SCALING_COMPILER=gcc \
      ./scripts/scaling_benchmark.sh inputFiles/benchmarkFiles/scaling/large.json 3
    '
  echo "Wrote results/scaling_linux.json"
  exit 0
fi

if command -v gh &>/dev/null; then
  echo "Triggering GitHub Actions workflow: scaling_linux.yml"
  gh workflow run scaling_linux.yml
  echo ""
  echo "When the run completes, download results with:"
  echo "  gh run list --workflow=scaling_linux.yml --limit 1"
  echo "  gh run download <run-id> -n scaling-linux-json -D results/"
  echo "  .venv/bin/python scripts/plot_scaling.py"
  exit 0
fi

echo "Install GitHub CLI (gh) or Docker to run Linux scaling." >&2
echo "  gh workflow run scaling_linux.yml" >&2
echo "  ./scripts/run_linux_scaling.sh --docker" >&2
exit 1
