# Known Issues and Limitations

## Field solver grid size

`numGrid` must be a **power of two** because the bundled FFT implementation requires it. Non-power-of-two grids will fail silently or produce incorrect fields.

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
- Thermal velocity initialization uses a non-seeded random number generator, so runs with `thermalVelocity > 0` are not reproducible.
- The UI assumes exactly two species with opposite drift velocities.

## Legacy input files

`inputFiles/landau.txt` and `inputFiles/twoStream.txt` use an old comma-separated format. The CLI accepts JSON only; use the JSON files under `inputFiles/` and `inputFiles/validation/` instead.

## Data type precision

`plasmaFrequency` and `chargeMassRatio` are stored as `double` in `DataStructs.h` but JSON inputs historically used integer values. Fractional values are supported.

## MPI

`CMakeLists.txt` previously referenced MPI include paths but MPI is not used. Parallelization (OpenMP/MPI) is planned work — see project roadmap in prior discussions.

## Energy conservation

The two-stream validation test allows up to 15% total-energy drift over 500 steps. This is a regression guard, not a tight physical conservation tolerance. Finer grids and smaller timesteps improve conservation.
