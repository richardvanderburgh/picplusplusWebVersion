#!/usr/bin/env python3
"""Generate the cold-plasma wave validation plot from a PIC++ JSON run."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputFiles" / "validation" / "coldPlasmaWave.json"
DEFAULT_OUTPUT = REPO_ROOT / "docs" / "images" / "cold_plasma_wave.png"


def find_binary() -> Path:
    candidates = [
        REPO_ROOT / "build" / "bin" / "PIC++Main",
        REPO_ROOT / "build" / "bin" / "Release" / "PIC++Main",
        REPO_ROOT / "build" / "bin" / "Release" / "PIC++Main.exe",
        REPO_ROOT / "buildLinux" / "bin" / "PIC++Main",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not find PIC++Main. Build the project first with Conan/CMake."
    )


def run_simulation(binary: Path, input_file: Path, results_file: Path) -> dict:
    results_file.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(binary), str(input_file), str(results_file)],
        check=True,
        cwd=REPO_ROOT,
    )
    return json.loads(results_file.read_text())


def rolling_envelope(values: np.ndarray, window: int) -> np.ndarray:
    env = np.empty_like(values)
    half = window // 2
    for i in range(len(values)):
        lo = max(0, i - half)
        hi = min(len(values), i + half + 1)
        env[i] = values[lo:hi].max()
    return env


def plot_results(data: dict, output_file: Path, dt: float) -> None:
    ese = np.asarray(data["ese"], dtype=float)
    ke = np.asarray(data["ke"][0], dtype=float)
    steps = np.arange(len(ese))
    time = steps * dt
    total = ke + ese

    envelope_window = max(15, int(round(2 * np.pi / dt)))
    env = rolling_envelope(ese, envelope_window)

    fig, axes = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)

    axes[0].plot(time, ese, color="#2980b9", linewidth=1, label="Grid field energy")
    axes[0].plot(time, env, color="#1a5276", linewidth=2, label="Rolling peak envelope")
    axes[0].set_ylabel("Field energy")
    axes[0].set_title("Cold plasma wave (v_th = 0): undamped oscillation")
    axes[0].legend(fontsize=8)
    axes[0].grid(True, alpha=0.3)
    axes[0].text(
        0.02,
        0.97,
        "omega ~ omega_p  =>  period T_p = 2 pi\n"
        "Integrated E^2 oscillates at ~2 omega_p",
        transform=axes[0].transAxes,
        va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    axes[1].plot(time, total, color="#27ae60", linewidth=1.5, label="Total energy")
    axes[1].set_xlabel("Time")
    axes[1].set_ylabel("Kinetic + field")
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.3)

    early = float(np.mean(ese[50:150]))
    late = float(np.mean(ese[-100:]))
    fig.text(
        0.5,
        0.01,
        f"Late/early mean field energy: {late / early:.2f} (expect ~1 for undamped cold wave)",
        ha="center",
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=150)
    print(f"Wrote validation plot to {output_file}")
    print(f"  late/early field-energy ratio: {late / early:.3f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results", type=Path, default=REPO_ROOT / "build" / "cold_plasma_results.json")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--binary", type=Path, default=None)
    args = parser.parse_args()

    binary = args.binary or find_binary()
    config = json.loads(args.input.read_text())
    dt = float(config.get("timeStepSize", 0.1))
    data = run_simulation(binary, args.input, args.results)
    plot_results(data, args.output, dt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
