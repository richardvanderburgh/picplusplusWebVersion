#!/usr/bin/env python3
"""Plot OpenMP strong-scaling results produced by scripts/scaling_benchmark.sh.

Reads results/scaling.json and writes docs/images/openmp_scaling.png, showing
measured speedup versus thread count against the ideal (linear) speedup.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results" / "scaling.json"
OUT = ROOT / "docs" / "images" / "openmp_scaling.png"


def main() -> int:
    if not RESULTS.exists():
        print(f"Missing {RESULTS}. Run scripts/scaling_benchmark.sh first.", file=sys.stderr)
        return 1

    data = json.loads(RESULTS.read_text())
    threads = data["threads"]
    times = data["time_loop_us"]

    if not threads or not times:
        print("No data in scaling.json.", file=sys.stderr)
        return 1

    baseline = times[0]
    speedup = [baseline / t for t in times]
    efficiency = [s / n for s, n in zip(speedup, threads)]

    fig, (ax_speed, ax_eff) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax_speed.plot(threads, threads, "--", color="gray", label="Ideal (linear)")
    ax_speed.plot(threads, speedup, "o-", color="#1f77b4", label="Measured")
    ax_speed.set_xlabel("OpenMP threads")
    ax_speed.set_ylabel("Speedup vs. 1 thread")
    ax_speed.set_title("PIC++ time-loop strong scaling")
    ax_speed.grid(True, alpha=0.3)
    ax_speed.legend()

    ax_eff.plot(threads, [1.0] * len(threads), "--", color="gray", label="Ideal")
    ax_eff.plot(threads, efficiency, "s-", color="#d62728", label="Measured")
    ax_eff.set_xlabel("OpenMP threads")
    ax_eff.set_ylabel("Parallel efficiency")
    ax_eff.set_ylim(0, 1.1)
    ax_eff.set_title("Parallel efficiency")
    ax_eff.grid(True, alpha=0.3)
    ax_eff.legend()

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130)
    print(f"Wrote {OUT}")

    print("\nthreads  time(us)     speedup  efficiency")
    for n, t, s, e in zip(threads, times, speedup, efficiency):
        print(f"{n:>7}  {t:>10}  {s:>8.2f}  {e:>9.2f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
