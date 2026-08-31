#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# Prefer project venv, then pip user install, then system PATH.
if [[ -x "$ROOT/.venv/bin/conan" ]]; then
  export PATH="$ROOT/.venv/bin:$PATH"
elif [[ -d "$HOME/Library/Python/3.9/bin" ]]; then
  export PATH="$HOME/Library/Python/3.9/bin:$PATH"
fi

if ! command -v conan &>/dev/null; then
  echo "Conan not found. Install with: python3 -m pip install 'conan>=2.0' cmake" >&2
  exit 1
fi

if ! conan --version 2>&1 | grep -q "Conan version 2"; then
  echo "Conan 2.x is required. Found: $(conan --version 2>&1)" >&2
  echo "Create a project venv with: python3 -m venv .venv && .venv/bin/pip install 'conan>=2.0' cmake" >&2
  exit 1
fi

if [[ -n "${PROFILE:-}" ]]; then
  :
elif [[ "$(uname -s)" == "Darwin" ]]; then
  PROFILE="buildUtils/macos_clang_release"
elif [[ "${OS:-}" == "Windows_NT" ]]; then
  PROFILE="buildUtils/win_release"
else
  PROFILE="buildUtils/linux_gcc_release"
fi

if [[ "$(uname -s)" == "Darwin" ]]; then
  libomp_found=false
  for prefix in /opt/homebrew/opt/libomp /usr/local/opt/libomp; do
    if [[ -f "$prefix/include/omp.h" ]]; then
      libomp_found=true
      break
    fi
  done
  if [[ "$libomp_found" == false ]]; then
    echo "Note: Homebrew libomp not found. The macOS build will be serial unless you run:" >&2
    echo "      brew install libomp" >&2
    echo "" >&2
  fi

  # Homebrew Qt is keg-only; export CMAKE_PREFIX_PATH so CMake finds Widgets + Charts.
  if [[ -z "${CMAKE_PREFIX_PATH:-}" ]] && command -v brew &>/dev/null; then
    for formula in qt qt@6; do
      if brew --prefix "$formula" &>/dev/null; then
        qt_prefix="$(brew --prefix "$formula")"
        if [[ -d "$qt_prefix/lib/cmake/Qt6" ]]; then
          export CMAKE_PREFIX_PATH="${qt_prefix}${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
          echo "Using Homebrew Qt at $qt_prefix"
          break
        fi
      fi
    done
  fi
fi

echo "Using Conan profile: $PROFILE"
conan install . --profile="$PROFILE" --build=missing -of=build
conan build . --profile="$PROFILE" -of=build

echo ""
echo "Build complete."
echo "  Run:  ./build/bin/PIC++Main inputFiles/exampleInput.json"
echo "  Test: ctest --test-dir build --output-on-failure"

if [[ -x "$ROOT/build/bin/PIC++GUI" ]]; then
  echo "  GUI:  ./build/bin/PIC++GUI"
elif [[ -d "$ROOT/gui" ]]; then
  echo "" >&2
  echo "PIC++GUI was not built (Qt 6 Widgets + Charts not found)." >&2
  if [[ "$(uname -s)" == "Darwin" ]]; then
    echo "  brew install qt" >&2
    echo "  export CMAKE_PREFIX_PATH=\"\$(brew --prefix qt)\"" >&2
    echo "  rm -f build/CMakeCache.txt && ./scripts/build.sh" >&2
  else
    echo "  Install Qt 6 development packages, then reconfigure CMake." >&2
  fi
elif [[ ! -d "$ROOT/gui" ]]; then
  echo "" >&2
  echo "PIC++GUI is not in this branch. Check out a branch with gui/ (e.g. cursor/3d-simulation-8672)." >&2
fi

if [[ -x "$ROOT/scripts/verify_openmp.sh" ]]; then
  echo ""
  if "$ROOT/scripts/verify_openmp.sh"; then
    :
  else
    echo "" >&2
    echo "OpenMP verification failed. Parallel speedups require an OpenMP-enabled build." >&2
    exit 1
  fi
fi
