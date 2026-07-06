#!/usr/bin/env python3
"""Plot Landau damping of a standing electron plasma wave in PIC++.

Runs the warm validation case (with spatial field snapshots) and the cold
control, then plots:
  1. E(x) at several times — standing-wave antinodes shrink (Landau damping)
  2. k=1 Fourier amplitude |E_k(t)| — mode envelope decay
  3. Warm vs. cold |E_k(t)| comparison
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputFiles" / "validation" / "landauDamping.json"
DEFAULT_COLD_INPUT = REPO_ROOT / "inputFiles" / "validation" / "coldPlasmaWave.json"
DEFAULT_OUTPUT = REPO_ROOT / "docs" / "images" / "landau_damping.png"


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


def grid_positions(num_points: int, length: float) -> np.ndarray:
    return np.linspace(0.0, length, num_points, endpoint=False)


def mode_amplitude(field: np.ndarray, x: np.ndarray, mode: int, length: float) -> float:
    phase = np.exp(-1j * 2.0 * np.pi * mode * x / length)
    return float(np.abs(np.sum(field * phase)) / field.size)


def mode_history(frames: list[dict], length: float, mode: int, dt: float, frame_period: int) -> tuple[np.ndarray, np.ndarray]:
    times = []
    amps = []
    for frame in frames:
        field = np.asarray(frame["electricField"], dtype=float)
        x = grid_positions(field.size, length)
        times.append(frame["frameNumber"] * dt)
        amps.append(mode_amplitude(field, x, mode, length))
    return np.asarray(times), np.asarray(amps)


def plot_results(
    warm: dict,
    cold: dict,
    output_file: Path,
    length: float,
    mode: int,
    dt: float,
    frame_period: int,
) -> None:
    warm_frames = warm.get("phaseFrames", [])
    if not warm_frames:
        raise RuntimeError("Warm run has no phaseFrames; set framePeriod > 0 in landauDamping.json")

    x = grid_positions(len(warm_frames[0]["electricField"]), length)
    warm_t, warm_k1 = mode_history(warm_frames, length, mode, dt, frame_period)

    cold_frames = cold.get("phaseFrames", [])
    if cold_frames:
        cold_t, cold_k1 = mode_history(cold_frames, length, mode, dt, frame_period)
    else:
        # Cold JSON uses framePeriod=0; approximate from integrated field energy envelope.
        cold_ese = np.asarray(cold["ese"], dtype=float)
        cold_t = np.arange(len(cold_ese)) * dt
        cold_k1 = np.sqrt(np.maximum(cold_ese, 1e-30))

    snapshot_indices = [0, len(warm_frames) // 4, len(warm_frames) // 2, 3 * len(warm_frames) // 4, -1]
    snapshot_labels = [f"t = {warm_frames[i]['frameNumber'] * dt:.1f}" for i in snapshot_indices]

    fig, axes = plt.subplots(3, 1, figsize=(9, 10))

    for idx, label in zip(snapshot_indices, snapshot_labels):
        field = np.asarray(warm_frames[idx]["electricField"], dtype=float)
        axes[0].plot(x, field, linewidth=1.2, alpha=0.85, label=label)
    axes[0].set_ylabel("E(x)")
    axes[0].set_title("Standing wave dissipates: field profiles at successive times (warm plasma)")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[0].grid(True, alpha=0.3)
    axes[0].text(
        0.02,
        0.97,
        "Nodes stay fixed; antinode amplitude shrinks\n"
        r"$E(x,t) \approx A(t)\cos(kx)\cos(\omega t)$",
        transform=axes[0].transAxes,
        va="top",
        fontsize=9,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    axes[1].plot(warm_t, warm_k1, "o-", color="#8e44ad", linewidth=1.5, markersize=3, label=r"$|E_{k=1}(t)|$ warm")
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"Mode amplitude $|E_{k=1}|$")
    axes[1].set_title("Landau damping: decaying k=1 envelope")
    axes[1].grid(True, alpha=0.3, which="both")
    axes[1].legend(fontsize=8)

    axes[2].plot(warm_t, warm_k1, "o-", color="#8e44ad", linewidth=1.5, markersize=3, label="Warm (v_th = 0.4)")
    axes[2].plot(cold_t, cold_k1, "o-", color="#2980b9", linewidth=1.5, markersize=3, label="Cold (v_th = 0)")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("Time")
    axes[2].set_ylabel(r"Mode amplitude")
    axes[2].set_title("Warm vs. cold control")
    axes[2].legend(fontsize=8)
    axes[2].grid(True, alpha=0.3, which="both")

    warm_ratio = warm_k1[-1] / warm_k1[0] if warm_k1[0] > 0 else 0.0
    cold_ratio = cold_k1[-1] / cold_k1[0] if cold_k1[0] > 0 else 0.0
    fig.text(
        0.5,
        0.01,
        f"|E_k| end/start — warm: {warm_ratio:.2f}, cold: {cold_ratio:.2f}. "
        "Landau damping: only the warm standing wave decays.",
        ha="center",
        fontsize=8,
    )

    fig.tight_layout(rect=(0, 0.03, 1, 1))
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=150)
    print(f"Wrote validation plot to {output_file}")
    print(f"  warm |E_k| end/start: {warm_ratio:.3f}")
    print(f"  cold |E_k| end/start: {cold_ratio:.3f}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--cold-input", type=Path, default=DEFAULT_COLD_INPUT)
    parser.add_argument("--results-warm", type=Path, default=REPO_ROOT / "build" / "landau_warm_results.json")
    parser.add_argument("--results-cold", type=Path, default=REPO_ROOT / "build" / "landau_cold_results.json")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--binary", type=Path, default=None)
    args = parser.parse_args()

    binary = args.binary or find_binary()
    config = json.loads(args.input.read_text())
    dt = float(config.get("timeStepSize", 0.1))
    length = float(config.get("spatialLength", 2 * np.pi))
    mode = int(config["species"][0].get("spatialPerturbationMode", 1))
    frame_period = int(config.get("framePeriod", 5))

    warm = run_simulation(binary, args.input, args.results_warm)

    cold_config = json.loads(args.cold_input.read_text())
    cold_config["framePeriod"] = frame_period
    cold_input = args.results_cold.with_suffix(".input.json")
    cold_input.write_text(json.dumps(cold_config, indent=2))
    cold = run_simulation(binary, cold_input, args.results_cold)

    plot_results(warm, cold, args.output, length, mode, dt, frame_period)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
