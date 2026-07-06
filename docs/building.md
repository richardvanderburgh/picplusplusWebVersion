# Building and Running PIC++

PIC++ uses [Conan 2.x](https://conan.io/) for dependencies and CMake for the build. The same Conan **profile** must be passed to both `conan install` and `conan build`; if they differ, dependency resolution or C++ standard settings can mismatch and the build will fail.

## Prerequisites

- C++20 compiler (GCC 10+, Clang 10+, Apple Clang, or MSVC 2019+)
- CMake 3.15+ (bundled with Conan or installed separately)
- Python 3.8+ with Conan 2.x

Recommended setup (avoids conflicts if Conan 1.x is also installed):

```bash
python3 -m venv .venv
.venv/bin/pip install "conan>=2.0" cmake
./scripts/build.sh   # uses .venv/bin/conan automatically
```

Or install globally:

```bash
python3 -m pip install "conan>=2.0" cmake
```

Ensure `conan` and `cmake` are on your `PATH`. `./scripts/build.sh` prepends `.venv/bin` or the common macOS pip user path automatically.

Verify you have Conan **2.x** (not 1.x):

```bash
conan --version   # should print "Conan version 2.x"
```

If you see a lockfile error, remove the legacy Conan 1 lockfile: `rm -f conan.lock`

## Quick build (recommended)

From the repo root:

```bash
./scripts/build.sh
```

This auto-selects a profile based on your OS:

| OS | Default profile |
|----|-----------------|
| macOS | `buildUtils/macos_clang_release` |
| Linux | `buildUtils/linux_gcc_release` |
| Windows | `buildUtils/win_release` |

Override with:

```bash
PROFILE=buildUtils/linux_clang_debug ./scripts/build.sh
```

## Manual build

### 1. Install dependencies

Pick the profile for your platform (all live under `buildUtils/`):

```bash
# Linux (GCC)
conan install . --profile=buildUtils/linux_gcc_release --build=missing -of=build

# macOS (Apple Clang)
conan install . --profile=buildUtils/macos_clang_release --build=missing -of=build

# Windows (MSVC)
conan install . --profile=buildUtils/win_release --build=missing -of=build
```

On macOS, if the profile's `compiler.version` does not match your Xcode version, either edit `buildUtils/macos_clang_release` or run `conan profile detect` and copy the detected values into that file.

### 2. Build and test

Pass the **same profile** to `conan build`:

```bash
conan build . --profile=buildUtils/macos_clang_release -of=build
```

`conan build` runs CTest automatically. To run tests again:

```bash
ctest --test-dir build --output-on-failure
```

### 3. Run a simulation

```bash
./build/bin/PIC++Main inputFiles/exampleInput.json
```

Write energies and phase frames to JSON:

```bash
./build/bin/PIC++Main inputFiles/validation/twoStreamInstability.json build/validation_results.json
```

On Windows the binary is at `build\bin\PIC++Main.exe` (or `build\bin\Release\PIC++Main.exe` depending on generator).

## Binaries

| Binary | Purpose |
|--------|---------|
| `build/bin/PIC++Main` | Run simulations from JSON input |
| `build/bin/PIC++Main_Test` | Integration, regression, and validation tests |
| `build/bin/PIC++Lib_Test` | Unit tests for library kernels |

Run a subset of tests:

```bash
./build/bin/PIC++Main_Test --gtest_filter="ValidationTest.*"
```

## Debug builds

```bash
conan install . --profile=buildUtils/linux_gcc_debug --build=missing -of=build
conan build . --profile=buildUtils/linux_gcc_debug -of=build
```

Profiles available: `linux_gcc_debug`, `linux_clang_debug`, `linux_clang_release`, `win_debug`, `win_release`, `macos_clang_release`.

## OpenMP (shared-memory parallelism)

The hot particle kernels are OpenMP-parallelized. OpenMP is **optional**: when
`find_package(OpenMP)` succeeds the parallel path is compiled in; otherwise the
code builds and runs correctly in serial. Linux GCC/Clang and MSVC ship OpenMP,
so `./scripts/build.sh` and CI pick it up automatically.

On macOS, install Homebrew's keg-only `libomp` once; `./scripts/build.sh` detects
it and configures CMake automatically:

```bash
brew install libomp
./scripts/build.sh
```

After every build, `scripts/build.sh` runs `scripts/verify_openmp.sh`, which
fails fast if the binary reports a serial build. CI runs the same check on
Ubuntu and Windows.

`PIC++Main` prints the active thread count at startup; control it with
`OMP_NUM_THREADS`. See [performance.md](performance.md) for the strong-scaling
study.

```bash
./scripts/scaling_benchmark.sh
.venv/bin/python scripts/plot_scaling.py
```

## Performance benchmarks

After building, run hyperfine sweeps (requires [hyperfine](https://github.com/sharkdp/hyperfine)):

```bash
./run_hyperfine.sh                          # timestep DOE (default)
./run_hyperfine.sh inputFiles/benchmarkFiles/particleDOE
```

Results are written to `results/`. Plot with `hyperfinePlottingTimestep.ipynb` or `hyperfineParticlePlotting.ipynb`.

## Validation plot

```bash
python3 -m pip install -r scripts/requirements.txt
python3 scripts/plot_two_stream_validation.py
```

See [validation.md](validation.md) for details.

## CI

GitHub Actions builds on Ubuntu (Docker) and Windows using `linux_gcc_release` and `win_release`, then runs `ctest --test-dir build --output-on-failure`.
