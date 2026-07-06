# Known Issues and Limitations

## Field solver grid size

`numGrid` must be a **power of two** because the bundled FFT implementation requires it.
Invalid values are rejected at startup with a clear error message from
`validateSimulationParams()` (CLI and `PICPlusPlus::initialize()`).

## Windows Conan profile / GitHub Actions runner drift

`buildUtils/win_release` and `buildUtils/win_debug` pin `compiler.version` to match whatever MSVC toolset is actually installed on the GitHub Actions `windows-latest` runner. GitHub periodically bumps the underlying Visual Studio version on that image (e.g. VS2022 -> VS2026), which changes the CMake generator Conan selects. If Windows CI fails with `CMake Error ... could not find any instance of Visual Studio`, run `conan profile detect` on a fresh `windows-latest` runner (or check the failing job's log for the `detect_api: Found msvc ...` / `Detected profile` output) and update `compiler.version` in both profiles to match.

## `conan build` requires `--profile` too

`conan build . -of=build` recomputes the dependency graph from scratch and, unlike `conan install`, does **not** remember the `--profile` passed to the preceding `conan install` call. Omitting `--profile` on the `conan build` step makes it silently fall back to Conan's auto-detected default profile, which can have a different `compiler.cppstd` (e.g. MSVC defaults to `14`) than the pinned build profiles (which use `20`). This produces a package-ID mismatch and an `ERROR: Missing binary` failure even though the preceding `conan install` succeeded. Always pass the same `--profile=...` to both `conan install` and `conan build` (see `scripts/build.sh` for the canonical pattern).

## `ctest` on Windows needs `-C <config>`

The Windows build uses the Visual Studio (multi-configuration) CMake generator, so the build type is chosen at test time rather than configure time. Running `ctest --test-dir build` without a config on Windows fails every test with `Test not available without configuration. (Missing "-C <config>"?)`. Pass `-C Release` (CI uses `-C ${{ env.BUILD_TYPE }}`). Linux/macOS use single-config generators (Make/Ninja) and don't need this.

## macOS Conan profile

`buildUtils/macos_clang_release` pins `compiler.version=14` and `arch=armv8` (Conan clamps Apple Clang to a max known version). Adjust these if you are on Intel macOS or a different Xcode/Clang version, or run `conan profile detect` and copy the detected values into the profile.

## Python environment (macOS `python3` vs `pip3` mismatch)

On macOS, bare `python3` often resolves to Homebrew's Python (e.g. `/opt/homebrew/bin/python3`, 3.14) while `pip3` targets the Command Line Tools Python (3.9). Installing a package with `pip3` can therefore land in a Python that `python3` doesn't use, producing `ModuleNotFoundError` even after a "successful" install.

Avoid this by using the project virtual environment for all Python work:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt        # Django UI
.venv/bin/python -m pip install -r scripts/requirements.txt # plotting
.venv/bin/python scripts/plot_two_stream_validation.py
.venv/bin/python manage.py runserver
```

The `.venv` also holds Conan/CMake, so a single environment drives builds, tests, plotting, and the web UI.

## Django UI

The web UI is a development aid, not a production front end:

- Hard-coded `ng = 32` in `templates/index.html` for plot axis limits does not track the form's `numGrid` value.
- Thermal velocity initialization uses a fixed seed (`std::mt19937(42)`), so runs
  with `thermalVelocity > 0` are reproducible across platforms.
- The UI assumes exactly two species with opposite drift velocities.

## Legacy input files

`inputFiles/landau.txt` and `inputFiles/twoStream.txt` use an old comma-separated format. The CLI accepts JSON only; use the JSON files under `inputFiles/` and `inputFiles/validation/` instead.

## Data type precision

`plasmaFrequency` and `chargeMassRatio` are stored as `double` in `DataStructs.h` but JSON inputs historically used integer values. Fractional values are supported.

## Parallelism (OpenMP / MPI)

The particle push, charge deposition, and kinetic-energy diagnostic are
OpenMP-parallelized (see [performance.md](performance.md)). Notes:

- **Floating-point determinism**: with OpenMP enabled, the charge-deposition
  merge and the energy reduction sum partial results in a thread-dependent
  order. Because floating-point addition is not associative, results can differ
  from the serial build (and between thread counts) at the ~1e-13 level. This is
  well within the regression/validation test tolerances, but bit-exact
  reproducibility across thread counts is not guaranteed.
- **MPI is not used.** Distributed-memory domain decomposition is future work;
  the stale MPI include paths were removed from `CMakeLists.txt`.

## Homebrew libomp is keg-only (macOS)

`brew install libomp` does not symlink into `/opt/homebrew`, so stock
`find_package(OpenMP)` fails for AppleClang. CMake now probes
`/opt/homebrew/opt/libomp` and `/usr/local/opt/libomp` automatically when you
run `./scripts/build.sh`. Without `libomp` installed the macOS build falls back
to serial (still correct, just single-threaded).

## Energy conservation

The two-stream validation test allows up to 15% total-energy drift over 500 steps. This is a regression guard, not a tight physical conservation tolerance. Finer grids and smaller timesteps improve conservation.
