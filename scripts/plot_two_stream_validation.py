#!/usr/bin/env python3
"""Generate the two-stream instability validation plot from a PIC++ JSON run."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputFiles" / "validation" / "twoStreamInstability.json"
DEFAULT_OUTPUT = REPO_ROOT / "docs" / "images" / "two_stream_validation.png"


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


def run_simulation(binary: Path, input_file: Path, results_file: Path) -> None:
    results_file.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(binary), str(input_file), str(results_file)],
        check=True,
        cwd=REPO_ROOT,
    )


def plot_results(results_file: Path, output_file: Path) -> None:
    with results_file.open() as handle:
        data = json.load(handle)

    ese = data["ese"]
    ke = data["ke"]
    time_steps = list(range(len(ese)))
    total_ke = [sum(species[step] for species in ke) for step in time_steps]
    total_energy = [ke_value + ese[step] for step, ke_value in enumerate(total_ke)]

    fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

    axes[0].plot(time_steps, ese, label="Electrostatic energy", color="#c0392b")
    axes[0].plot(time_steps, total_ke, label="Total kinetic energy", color="#2980b9")
    axes[0].set_ylabel("Energy")
    axes[0].set_title("Two-Stream Instability Validation")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(time_steps, total_energy, label="Total energy", color="#27ae60")
    axes[1].set_xlabel("Time step")
    axes[1].set_ylabel("Total energy")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    fig.tight_layout()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=150)
    print(f"Wrote validation plot to {output_file}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--results", type=Path, default=REPO_ROOT / "build" / "validation_results.json")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--binary", type=Path, default=None)
    args = parser.parse_args()

    binary = args.binary or find_binary()
    run_simulation(binary, args.input, args.results)
    plot_results(args.results, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
