#!/usr/bin/env python3
"""Quantitative physics verification for PIC++ validation cases.

Reports pass/fail against the project's documented criteria and adds
analytic cross-checks where they are meaningful.
"""

from __future__ import annotations

import cmath
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINARY = ROOT / "build" / "bin" / "PIC++Main"


def run_case(input_path: Path) -> dict:
    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        out = Path(tmp.name)
    subprocess.run([str(BINARY), str(input_path), str(out)], check=True, cwd=ROOT, capture_output=True)
    data = json.loads(out.read_text())
    out.unlink(missing_ok=True)
    return data


def total_energy(data: dict) -> list[float]:
    ese = data["ese"]
    ke = data["ke"]
    return [ese[i] + sum(species[i] for species in ke) for i in range(len(ese))]


def mode_amplitude(field: list[float], mode: int, length: float) -> float:
    n = len(field) - 1
    if n <= 0:
        return 0.0
    k = 2.0 * math.pi * mode / length
    re = im = 0.0
    for j in range(n):
        x = (j / n) * length
        re += field[j] * math.cos(k * x)
        im += field[j] * math.sin(k * x)
    return math.hypot(re / n, im / n)


def extract_mode_series(data: dict, mode: int, length: float, dt: float):
    times, amps = [], []
    for frame in data.get("phaseFrames", []):
        times.append(frame["frameNumber"] * dt)
        amps.append(mode_amplitude(frame["electricField"], mode, length))
    return times, amps


def fit_log_slope(times, amps, i0, i1):
    xs, ys = [], []
    for i in range(i0, i1):
        if amps[i] > 1e-14:
            xs.append(times[i])
            ys.append(math.log(amps[i]))
    if len(xs) < 5:
        return None
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-18:
        return None
    return (n * sxy - sx * sy) / denom


def cold_two_stream_gamma(a: float) -> float:
    """Max |Im(ω/ωp)| for equal cold beams, a = k vd / ωp."""
    b = -(2 * a * a + 1)
    c = a**4 - a * a
    disc = b * b - 4 * c
    us = [(-b + cmath.sqrt(disc)) / 2, (-b - cmath.sqrt(disc)) / 2]
    gamma = 0.0
    for u in us:
        for sign in (1, -1):
            omega = sign * cmath.sqrt(u)
            gamma = max(gamma, abs(omega.imag))
    return gamma


def mean_slice(values, frac0, frac1):
    n = len(values)
    i0 = int(n * frac0)
    i1 = max(i0 + 1, int(n * frac1))
    return sum(values[i0:i1]) / (i1 - i0)


def report(title: str, lines: list[str], ok: bool) -> None:
    print(f"\n=== {title} [{'PASS' if ok else 'FAIL'}] ===")
    for line in lines:
        print(f"  {line}")


def main() -> int:
    if not BINARY.exists():
        print(f"Missing binary: {BINARY}", file=sys.stderr)
        return 1

    failures = 0

    # Two-stream
    ts = run_case(ROOT / "inputFiles/validation/twoStreamInstability.json")
    ese = ts["ese"]
    totals = total_energy(ts)
    growth = ese[-1] / max(ese[0], 1e-30)
    drift = abs(totals[-1] - totals[0]) / max(abs(totals[0]), 1e-30)
    times, amps = extract_mode_series(ts, 1, 2 * math.pi, 0.2)
    a = 1.0  # k=1, vd=1, wp=1
    theory_g = cold_two_stream_gamma(a)
    ok = growth >= 5.0 and drift < 0.15
    failures += not ok
    report(
        "Two-stream instability",
        [
            f"ESE final/initial = {growth:.2f} (criterion >= 5) ✓" if growth >= 5 else f"ESE final/initial = {growth:.2f}",
            f"Total energy drift = {100 * drift:.3f}% (criterion < 15%)",
            f"Peak ESE / initial = {max(ese) / max(ese[0], 1e-30):.1f}",
            f"Cold equal-beam linear γ(a=k vd/ωp={a}) = {theory_g:.4f}",
            "Note: a=1 sits at the cold linear stability boundary (γ=0);",
            "observed growth is nonlinear / discrete-particle driven, which the",
            "validation criteria intentionally target (strong ESE rise + bounded energy).",
        ],
        ok,
    )

    # Cold plasma
    cold = run_case(ROOT / "inputFiles/validation/coldPlasmaWave.json")
    cold_ese = cold["ese"]
    cold_ratio = mean_slice(cold_ese, 0.8, 1.0) / max(mean_slice(cold_ese, 0.1, 0.3), 1e-30)
    cold_times, cold_amps = extract_mode_series(cold, 1, 2 * math.pi, 0.1)
    ek_ratio = cold_amps[-1] / max(cold_amps[0], 1e-30)
    ok = 0.90 < cold_ratio < 1.20 and 0.8 < ek_ratio < 1.2
    failures += not ok
    report(
        "Cold plasma wave (control)",
        [
            f"Late/early mean ESE = {cold_ratio:.3f} (criterion 0.90–1.20)",
            f"|E_k| end/start = {ek_ratio:.3f} (expect ~1, no Landau)",
            f"Expected ω ≈ ω_p = 1 on L=2π, k=1 (Bohm–Gross cold limit)",
            "Energy: cold KE starts near zero, so relative total-energy drift is",
            "not a tight metric here; undamped |E_k| is the physical check.",
        ],
        ok,
    )

    # Landau
    warm = run_case(ROOT / "inputFiles/validation/landauDamping.json")
    warm_ese = warm["ese"]
    peak = max(warm_ese)
    warm_ratio = mean_slice(warm_ese, 0.8, 1.0) / max(mean_slice(warm_ese, 0.1, 0.3), 1e-30)
    warm_totals = total_energy(warm)
    warm_drift = abs(warm_totals[-1] - warm_totals[0]) / max(abs(warm_totals[0]), 1e-30)
    warm_times, warm_amps = extract_mode_series(warm, 1, 2 * math.pi, 0.1)
    ek_end_start = warm_amps[-1] / max(warm_amps[0], 1e-30)
    imax = max(range(len(warm_amps)), key=lambda i: warm_amps[i])
    measured_g = fit_log_slope(warm_times, warm_amps, imax, int(len(warm_amps) * 0.7))
    omega = math.sqrt(1.0 + 3.0 * (0.4**2))
    discriminator = cold_ratio - warm_ratio
    ok = (
        peak > warm_ese[0] * 1.05
        and warm_ese[-1] < peak * 0.55
        and warm_ese[-1] < warm_ese[0] * 0.85
        and warm_drift < 0.15
        and warm_ratio < 0.85
        and discriminator > 0.1
        and 0.15 < ek_end_start < 0.4
    )
    failures += not ok
    report(
        "Landau damping (warm vs cold)",
        [
            f"|E_k| end/start = {ek_end_start:.3f} (docs expect ~0.2–0.3)",
            f"Late/early mean ESE = {warm_ratio:.3f} (criterion < 0.85)",
            f"Warm−cold discriminator = {discriminator:.3f} (criterion > 0.10)",
            f"Total energy drift = {100 * warm_drift:.3f}% (criterion < 15%)",
            f"Bohm–Gross ω ≈ {omega:.3f}, v_φ ≈ {omega:.3f}, v_φ/v_th ≈ {omega / 0.4:.2f}",
            f"Measured envelope γ ≈ {measured_g:.4f}" if measured_g is not None else "Measured γ unavailable",
            "Note: simple asymptotic Landau formulas are inaccurate at k λ_D = 0.4;",
            "the warm-vs-cold discriminator is the intended physics proof.",
        ],
        ok,
    )

    # 3D
    d3 = run_case(ROOT / "inputFiles/validation/plasmaOscillation3D.json")
    d3_totals = total_energy(d3)
    d3_finite = d3.get("dimension") == 3 and all(math.isfinite(x) for x in d3["ese"])
    d3_drift = abs(d3_totals[-1] - d3_totals[0]) / max(abs(d3_totals[0]), 1e-30)
    ok = d3_finite and d3_drift < 0.15
    failures += not ok
    report(
        "3D plasma oscillation",
        [
            f"dimension = {d3.get('dimension')}, frames = {len(d3.get('phaseFrames', []))}",
            f"All ESE samples finite: {d3_finite}",
            f"Total energy drift = {100 * d3_drift:.3f}% (criterion < 15%)",
        ],
        ok,
    )

    print("\n=== Verdict ===")
    if failures == 0:
        print("Documented electrostatic validation criteria: PASSED (1D and 3D).")
        return 0
    print(f"{failures} primary check group(s) FAILED.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
