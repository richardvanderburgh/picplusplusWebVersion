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

Four kernels run every time step. The particle kernels scale with particle count;
the field solve scales with grid size. Each is parallelized according to its
data-dependency pattern.

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

### 4. Field solve (`fields`) — k-space Poisson and grid derivatives

After the charge density FFT (serial butterfly in `CFFT`), the Poisson solve in
Fourier space and the electric-field finite-difference stencils are independent
per mode/grid point and run as parallel `for` loops:

```cpp
#pragma omp parallel for schedule(static)
for (int k = 0; k < numGrid + 1; ++k) {
    // solve for complexPotentialK[k]
}

#pragma omp parallel for schedule(static)
for (int j = 1; j < numGrid; ++j) {
    electricField[j] = (phi[j - 1] - phi[j + 1]) / (2 * dx);
}
```

The FFT itself remains serial for now; parallelizing its butterfly stages needs
cache-aware tiling to avoid false sharing on medium grids (512 cells).

All OpenMP pragmas are guarded with `#ifdef _OPENMP`, so a build without OpenMP
(e.g. stock AppleClang) compiles the identical code path serially and correctly
— the parallelization is purely opt-in at configure time.

## Strong-scaling results

Fixed problem size (**800,000 particles**, 512 grid cells, 200 steps;
`framePeriod: 0` so frame I/O doesn't dominate), varying the OpenMP thread count.
Measured on **macOS / AppleClang** (Apple Silicon, 14 logical cores) and
**Linux / GCC** (GitHub Actions `ubuntu-latest`, 4 vCPUs).

![OpenMP strong scaling](images/openmp_scaling.png)

### macOS / AppleClang (local)

| Threads | Time loop (µs) | Speedup | Efficiency |
|--------:|---------------:|--------:|-----------:|
| 1       | 642,165        | 1.00×   | 100%       |
| 2       | 339,097        | 1.89×   | 95%        |
| 4       | 206,914        | 3.10×   | 78%        |
| 8       | 147,396        | 4.36×   | 54%        |
| 14      | 163,271        | 3.93×   | 28%        |

### Linux / GCC (CI)

| Threads | Time loop (µs) | Speedup | Efficiency |
|--------:|---------------:|--------:|-----------:|
| 1       | 1,037,396      | 1.00×   | 100%       |
| 2       | 529,205        | 1.96×   | 98%        |
| 4       | 440,881        | 2.35×   | 59%        |

Raw JSON: `results/scaling.json` (local) and `results/scaling_linux.json` (CI).

Observations:

- Both platforms show near-linear speedup at 2 threads (~95–98% efficiency).
- macOS continues scaling to **4.36× at 8 threads**, peaking near the
  performance-core count; efficiency drops past 8 threads as work spills onto
  slower efficiency cores.
- Linux GCC reaches **2.35× at 4 threads** on the 4-vCPU CI runner. Absolute
  single-thread time is higher than Apple Silicon (expected for a shared
  virtualized runner), and efficiency at 4 threads is lower as the serial FFT
  fraction and memory bandwidth become the bottleneck on a smaller socket.
- The residual serial fraction — the FFT-based Poisson field solve
  (`fields`, O(numGrid log numGrid)) and per-step diagnostics — caps the
  achievable speedup per Amdahl's law. On a bandwidth-bound streaming kernel
  like this, sub-linear scaling is expected and honest.

### Next steps

- Parallelize the FFT butterfly stages (attempted; needs cache-aware tiling to
  beat serial on medium grids).
- NUMA-aware first-touch initialization on multi-socket nodes.
- Distributed-memory domain decomposition (MPI) for multi-node runs.
- SIMD intrinsics or explicit `#pragma omp simd` on the gather in `accel`.

## Reproducing the benchmark

An OpenMP-enabled build is required (Linux GCC/Clang ship it; on macOS install
`libomp` via `brew install libomp`).

```bash
# Build (see docs/building.md; on macOS: brew install libomp first).
./scripts/build.sh

# Run the strong-scaling sweep and plot it.
./scripts/scaling_benchmark.sh
.venv/bin/python scripts/plot_scaling.py
```

Results are written to `results/scaling.json` and the plot to
`docs/images/openmp_scaling.png`.

### Linux GCC benchmark (CI)

The workflow `.github/workflows/scaling_linux.yml` builds with GCC inside the
project's Ubuntu container image and uploads `results/scaling_linux.json`.
Trigger it manually or let it run on pushes to `main` that touch the solver:

```bash
gh workflow run scaling_linux.yml
gh run list --workflow=scaling_linux.yml --limit 1
gh run download <run-id> -n scaling-linux-json -D results/
.venv/bin/python scripts/plot_scaling.py
```

Or use the helper: `./scripts/run_linux_scaling.sh` (wraps `gh workflow run`).

### Verifying OpenMP is enabled

`./scripts/build.sh` runs this automatically after configure/build:

```bash
scripts/verify_openmp.sh
```

CI runs the same check on Ubuntu and Windows. A passing run prints
`OK: OpenMP enabled (...)`; a serial build fails with install/rebuild guidance.
