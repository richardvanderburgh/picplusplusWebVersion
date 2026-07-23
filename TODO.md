# Roadmap

Living checklist of project status. Prefer this over open GitHub issues for
current intent; several older issues are stale relative to `main`.

## Done

- JSON input files and CLI (`inputFiles/`, validation, demos)
- GoogleTest unit, regression, and physics validation suites
- Django web UI with Plotly plots (phase space, E(x), energy, \|Eₖ\|)
- Demo catalog + parameter form + JSON/CSV file import
- OpenMP parallelization of particle push, charge deposition, and KE diagnostics
- CI on Ubuntu (Docker) and Windows; Linux OpenMP scaling workflow
- Landau damping / two-stream / cold-plasma validation cases and docs

## Near-term

- [ ] Close or update stale GitHub issues that are already implemented
  (notably #23 OpenMP; review #17 data export vs current JSON output)
- [ ] Per-species editors in the web form (today: shared template, or full JSON import)
- [ ] Optional configurable RNG seed for thermal velocity (default remains fixed for reproducibility)
- [ ] Tighten energy-conservation validation on finer grids / smaller Δt where useful

## Larger / exploratory

- [ ] Parallelize the FFT field solve (cache-aware tiling; see `docs/performance.md`)
- [ ] SIMD / `#pragma omp simd` on the gather in `accel`
- [ ] NUMA-aware first-touch initialization on multi-socket nodes
- [ ] Distributed-memory domain decomposition (MPI) — GitHub #22
- [ ] Optional FFTW3 backend instead of the bundled power-of-two FFT — GitHub #12
- [ ] Desktop rendering window class (SDL) — GitHub #19
- [ ] DOE automation shell scripts for source optimization — GitHub #25
- [ ] React front end (optional rewrite; current UI is Django + Plotly)
- [ ] Emscripten / in-browser WASM port (prior spike hit pointer issues on larger inputs)

## Not planned right now

- Production hardening of the Django UI (auth, multi-user jobs, etc.) — it remains a local development aid
