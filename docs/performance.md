# Performance & Parallelism

PIC++ is a 1D electrostatic particle-in-cell code. Its cost is dominated by the
per-particle work in the time-integration loop, so the performance strategy
focuses on (1) a cache-friendly data layout and (2) shared-memory
parallelization of the hot particle kernels with OpenMP.

## Data layout: Structure-of-Arrays (SoA)

The compute-critical state is stored as a **Structure-of-Arrays**. Each species
keeps its particle data in separate contiguous `std::vector<double>` arrays
rather than an array of particle structs:

```cpp
struct SpeciesData {
    // ...
    std::vector<double> particlePositions;    // x_0, x_1, x_2, ...
    std::vector<double> particleXVelocities;  // v_0, v_1, v_2, ...
};
```

Why this matters for a PIC push:

- The particle push and charge deposition stream linearly through positions and
  velocities. SoA keeps each stream contiguous, maximizing cache-line
  utilization and enabling the compiler to auto-vectorize the arithmetic.
- An Array-of-Structures layout (`struct Particle { pos; vel; species; id; }`)
  would interleave fields the hot loop doesn't use (`species`, `id`), wasting
  bandwidth on every cache line.

The AoS `Particle` struct still exists, but only for **diagnostic output**
(serializing phase-space frames to JSON) — never on the compute path.

## OpenMP parallelization of the hot loops

Three kernels run every time step and scale with particle count. Each is
parallelized according to its data-dependency pattern.

### 1. Particle push (`accel`) — embarrassingly parallel

Every particle reads the shared acceleration grid (read-only) and writes only
its own velocity. There are no cross-particle dependencies, so it is a plain
parallel `for`:

```cpp
#pragma omp parallel for schedule(static)
for (int i = 0; i < numParticles; ++i) {
    // gather grid acceleration, update velocities[i]
}
```

### 2. Charge deposition (`move`) — scatter with a write race

This is the classic PIC parallelization challenge. Each particle *scatters*
charge into two grid cells:

```cpp
rho[j]     += q - drho;
rho[j + 1] += drho;
```

Two particles in the same cell would race on `rho[j]`. PIC++ resolves this
**without atomics** by giving each thread a private density buffer, then merging
the buffers once at the end of the parallel region:

```cpp
#pragma omp parallel
{
    std::vector<double> localRho(rho.size(), 0.0);

    #pragma omp for nowait schedule(static)
    for (int i = 0; i < numParticles; ++i) {
        // update position (independent per particle)
        localRho[j]     += q - drho;   // no contention: thread-private
        localRho[j + 1] += drho;
    }

    #pragma omp critical
    for (size_t k = 0; k < rho.size(); ++k) rho[k] += localRho[k];  // merge
}
```

This trades a small amount of memory (one grid-sized buffer per thread; the grid
is tiny compared to the particle count) for a race-free deposition that avoids
the heavy atomic contention a naive `#pragma omp atomic` would incur when many
particles land in the same cell.

### 3. Kinetic-energy diagnostic — reduction

The kinetic energy is a sum over all particles, expressed as an OpenMP
reduction:

```cpp
double kineticEnergy = 0.0;
#pragma omp parallel for reduction(+ : kineticEnergy) schedule(static)
for (int i = 0; i < numParticles; ++i)
    kineticEnergy += 0.5 * v[i] * v[i] * mass;
```

All OpenMP pragmas are guarded with `#ifdef _OPENMP`, so a build without OpenMP
(e.g. stock AppleClang) compiles the identical code path serially and correctly
— the parallelization is purely opt-in at configure time.

## Strong-scaling results

Fixed problem size (**800,000 particles**, 512 grid cells, 200 steps;
`framePeriod: 0` so frame I/O doesn't dominate), varying the OpenMP thread count.
Measured on Apple Silicon (8 performance + efficiency cores).

![OpenMP strong scaling](images/openmp_scaling.png)

| Threads | Time loop (µs) | Speedup | Efficiency |
|--------:|---------------:|--------:|-----------:|
| 1       | 644,128        | 1.00×   | 100%       |
| 2       | 346,665        | 1.86×   | 93%        |
| 4       | 208,176        | 3.09×   | 77%        |
| 8       | 149,458        | 4.31×   | 54%        |
| 14      | 166,370        | 3.87×   | 28%        |

Observations:

- Near-linear speedup through 2–4 threads (77–93% efficiency).
- Speedup peaks around the physical performance-core count (~8) and then
  regresses as work spills onto slower efficiency cores and the run becomes
  memory-bandwidth bound.
- The residual serial fraction — the FFT-based Poisson field solve
  (`fields`, O(numGrid log numGrid)) and per-step diagnostics — caps the
  achievable speedup per Amdahl's law. On a bandwidth-bound streaming kernel
  like this, sub-linear scaling is expected and honest.

### Next steps

- Parallelize / vectorize the FFT field solve to shrink the serial fraction.
- NUMA-aware first-touch initialization on multi-socket nodes.
- Distributed-memory domain decomposition (MPI) for multi-node runs.
- SIMD intrinsics or explicit `#pragma omp simd` on the gather in `accel`.

## Reproducing the benchmark

An OpenMP-enabled build is required (Linux GCC/Clang ship it; on macOS install
`libomp` via `brew install libomp`).

```bash
# Build (see docs/building.md; on macOS pass the libomp hints so
# find_package(OpenMP) succeeds with AppleClang).
./scripts/build.sh

# Run the strong-scaling sweep and plot it.
PICPP_BIN=build/bin/PIC++Main ./scripts/scaling_benchmark.sh
.venv/bin/python scripts/plot_scaling.py
```

Results are written to `results/scaling.json` and the plot to
`docs/images/openmp_scaling.png`.

### macOS + AppleClang note

Homebrew's `libomp` is keg-only, so CMake needs explicit hints. Configure with:

```bash
cmake -S . -B build-omp \
  -DCMAKE_TOOLCHAIN_FILE="$PWD/build/conan_toolchain.cmake" \
  -DCMAKE_BUILD_TYPE=Release \
  -DOpenMP_ROOT=/opt/homebrew/opt/libomp \
  -DOpenMP_CXX_FLAGS="-Xclang -fopenmp -I/opt/homebrew/opt/libomp/include" \
  -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib
cmake --build build-omp -j
```

(Run `./scripts/build.sh` first once to generate `build/conan_toolchain.cmake`.)
