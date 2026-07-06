#!/usr/bin/env python3
"""Plot OpenMP strong-scaling results from scripts/scaling_benchmark.sh.

Reads results/scaling.json (local platform) and, when present,
results/scaling_linux.json (GCC/Linux CI or Docker run). Writes
docs/images/openmp_scaling.png.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RESULT_FILES = [
    ("Local", ROOT / "results" / "scaling.json"),
    ("Linux GCC", ROOT / "results" / "scaling_linux.json"),
]
OUT = ROOT / "docs" / "images" / "openmp_scaling.png"


def load_run(path: Path) -> dict | None:
    if not path.exists():
        return None
    data = json.loads(path.read_text())
    threads = data.get("threads", [])
    times = data.get("time_loop_us", [])
    if not threads or not times:
        return None
    platform = data.get("platform", "unknown")
    compiler = data.get("compiler", "unknown")
    nproc = data.get("nproc", "?")
    label = f"{platform} / {compiler} ({nproc} cores)"
    baseline = times[0]
    speedup = [baseline / t for t in times]
    efficiency = [s / n for s, n in zip(speedup, threads)]
    return {
        "label": label,
        "threads": threads,
        "times": times,
        "speedup": speedup,
        "efficiency": efficiency,
    }


def main() -> int:
    runs: list[tuple[str, dict]] = []
    for name, path in RESULT_FILES:
        loaded = load_run(path)
        if loaded is not None:
            runs.append((name, loaded))

    if not runs:
        print(
            "No scaling data found. Run scripts/scaling_benchmark.sh and/or "
            "scripts/run_linux_scaling.sh",
            file=sys.stderr,
        )
        return 1

    fig, (ax_speed, ax_eff) = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

    max_threads = max(max(r["threads"]) for _, r in runs)
    ax_speed.plot(
        [1, max_threads],
        [1, max_threads],
        "--",
        color="gray",
        linewidth=1.5,
        label="Ideal (linear)",
    )

    for idx, (_, run) in enumerate(runs):
        color = colors[idx % len(colors)]
        ax_speed.plot(
            run["threads"],
            run["speedup"],
            "o-",
            color=color,
            linewidth=2,
            markersize=6,
            label=run["label"],
        )
        ax_eff.plot(
            run["threads"],
            run["efficiency"],
            "s-",
            color=color,
            linewidth=2,
            markersize=6,
            label=run["label"],
        )

    ax_speed.set_xlabel("OpenMP threads")
    ax_speed.set_ylabel("Speedup vs. 1 thread")
    ax_speed.set_title("PIC++ time-loop strong scaling")
    ax_speed.grid(True, alpha=0.3)
    ax_speed.legend(fontsize=8)

    ax_eff.axhline(1.0, color="gray", linestyle="--", linewidth=1.5, label="Ideal")
    ax_eff.set_xlabel("OpenMP threads")
    ax_eff.set_ylabel("Parallel efficiency")
    ax_eff.set_ylim(0, 1.1)
    ax_eff.set_title("Parallel efficiency")
    ax_eff.grid(True, alpha=0.3)
    ax_eff.legend(fontsize=8)

    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=130)
    print(f"Wrote {OUT}")

    for name, run in runs:
        print(f"\n{name} ({run['label']})")
        print("threads  time(us)     speedup  efficiency")
        for n, t, s, e in zip(run["threads"], run["times"], run["speedup"], run["efficiency"]):
            print(f"{n:>7}  {t:>10}  {s:>8.2f}  {e:>9.2f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
