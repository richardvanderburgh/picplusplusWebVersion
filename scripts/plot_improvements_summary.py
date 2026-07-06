#!/usr/bin/env python3
"""Generate a multi-panel summary of PIC++ Tier 1/2 improvements.

Panels:
  1. OpenMP strong-scaling speedup
  2. Two-stream instability energy growth (physics validation)
  3. Wall-clock time vs thread count (absolute)
  4. Energy budget over the two-stream run
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
SCALING = ROOT / "results" / "scaling.json"
VALIDATION_INPUT = ROOT / "inputFiles" / "validation" / "twoStreamInstability.json"
VALIDATION_JSON = ROOT / "build" / "validation_results.json"
OUT = ROOT / "docs" / "images" / "improvements_summary.png"


def find_binary() -> Path:
    for candidate in (
        ROOT / "build" / "bin" / "PIC++Main",
        ROOT / "build" / "bin" / "Release" / "PIC++Main",
        ROOT / "build" / "bin" / "Release" / "PIC++Main.exe",
    ):
        if candidate.exists():
            return candidate
    raise FileNotFoundError("PIC++Main not found — run ./scripts/build.sh first")


def ensure_validation_data(binary: Path) -> dict:
    VALIDATION_JSON.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(binary), str(VALIDATION_INPUT), str(VALIDATION_JSON)],
        check=True,
        cwd=ROOT,
    )
    return json.loads(VALIDATION_JSON.read_text())


def main() -> int:
    if not SCALING.exists():
        print(f"Missing {SCALING}. Run scripts/scaling_benchmark.sh first.", file=sys.stderr)
        return 1

    scaling = json.loads(SCALING.read_text())
    threads = scaling["threads"]
    times_us = scaling["time_loop_us"]
    times_s = [t / 1e6 for t in times_us]
    baseline = times_us[0]
    speedup = [baseline / t for t in times_us]

    binary = find_binary()
    validation = ensure_validation_data(binary)
    ese = validation["ese"]
    ke = validation["ke"]
    steps = list(range(len(ese)))
    total_ke = [sum(species[step] for species in ke) for step in steps]
    total_energy = [k + ese[step] for step, k in enumerate(total_ke)]

  # 2x2 dashboard
    fig = plt.figure(figsize=(12, 9))
    fig.suptitle("PIC++ improvements summary", fontsize=14, fontweight="bold", y=0.98)

    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(threads, threads, "--", color="#999", linewidth=1.5, label="Ideal")
    ax1.plot(threads, speedup, "o-", color="#1f77b4", linewidth=2, markersize=7, label="Measured")
    ax1.set_xlabel("OpenMP threads")
    ax1.set_ylabel("Speedup vs. 1 thread")
    ax1.set_title("Strong scaling (800k particles)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper left")
    peak_idx = max(range(len(speedup)), key=lambda i: speedup[i])
    ax1.annotate(
        f"peak {speedup[peak_idx]:.2f}× @ {threads[peak_idx]} threads",
        xy=(threads[peak_idx], speedup[peak_idx]),
        xytext=(threads[peak_idx] + 1, speedup[peak_idx] - 0.6),
        arrowprops=dict(arrowstyle="->", color="#333"),
        fontsize=9,
    )

    ax2 = fig.add_subplot(2, 2, 2)
    bars = ax2.bar(
        [str(t) for t in threads],
        times_s,
        color=["#aec7e8" if t == 1 else "#1f77b4" for t in threads],
        edgecolor="#333",
    )
    ax2.set_xlabel("OpenMP threads")
    ax2.set_ylabel("Time loop (seconds)")
    ax2.set_title("Absolute runtime (lower is better)")
    ax2.grid(True, axis="y", alpha=0.3)
    for bar, t_s in zip(bars, times_s):
        ax2.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.02,
            f"{t_s:.2f}s",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax3 = fig.add_subplot(2, 2, 3)
    ax3.plot(steps, ese, color="#c0392b", linewidth=1.5, label="Field energy")
    ax3.plot(steps, total_ke, color="#2980b9", linewidth=1.5, label="Kinetic energy")
    ax3.set_xlabel("Time step")
    ax3.set_ylabel("Energy")
    ax3.set_title("Two-stream instability (physics validation)")
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    growth = ese[-1] / max(ese[0], 1e-12)
    ax3.text(
        0.03,
        0.95,
        f"Field energy growth: {growth:.1f}×",
        transform=ax3.transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    ax4 = fig.add_subplot(2, 2, 4)
    drift = abs(total_energy[-1] - total_energy[0]) / total_energy[0] * 100
    ax4.plot(steps, total_energy, color="#27ae60", linewidth=1.5)
    ax4.set_xlabel("Time step")
    ax4.set_ylabel("Total energy")
    ax4.set_title("Energy budget conservation")
    ax4.grid(True, alpha=0.3)
    ax4.axhline(total_energy[0], color="#999", linestyle="--", linewidth=1, label="Initial")
    ax4.text(
        0.03,
        0.95,
        f"Total drift: {drift:.1f}%",
        transform=ax4.transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )
    ax4.legend()

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=150)
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
