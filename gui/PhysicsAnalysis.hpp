#pragma once

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include "SimulationConfig.h"

namespace PhysicsAnalysis {

struct GrowthFit {
	double rate = 0.0;
	double amplitude0 = 0.0;
	double tStart = 0.0;
	double tEnd = 0.0;
	bool valid = false;
	std::string label;
};

struct TheoryCurve {
	std::string name;
	double rate = 0.0; // γ for exp(γ t) or -γ for damping
	bool isDamping = false;
	bool valid = false;
	std::string note;
};

inline double waveNumber(const FormParams& params) {
	if (params.spatialLength <= 0.0) {
		return 0.0;
	}
	return (2.0 * 3.14159265358979323846 * params.spatialPerturbationMode) / params.spatialLength;
}

// Approximate Landau damping rate for a warm Maxwellian (fluid-ish estimate).
// γ ≈ −√(π/8) · (ωp³ / (k³ vth³)) · exp(−1/(2 k² λd²) − 3/2) with λd = vth/ωp.
// Used as a teaching overlay, not a high-fidelity kinetic solver.
inline TheoryCurve landauDampingEstimate(const FormParams& params) {
	TheoryCurve curve;
	curve.name = "Landau estimate";
	curve.isDamping = true;
	const double k = waveNumber(params);
	const double vth = params.thermalVelocity;
	const double wp = params.plasmaFrequency;
	if (k <= 0.0 || vth <= 1e-8 || wp <= 0.0) {
		curve.note = "Need finite v_th and mode k > 0 for Landau estimate.";
		return curve;
	}
	const double kLd = k * vth / wp;
	if (kLd < 0.05) {
		curve.note = "k λ_D too small — damping estimate unreliable.";
		return curve;
	}
	const double expArg = -0.5 / (kLd * kLd) - 1.5;
	const double gamma = -std::sqrt(3.14159265358979323846 / 8.0) * (wp / (kLd * kLd * kLd)) * std::exp(expArg);
	curve.rate = gamma; // negative
	curve.valid = std::isfinite(gamma) && gamma < 0.0;
	curve.note = "Approx γ for Maxwellian Landau damping (teaching overlay).";
	return curve;
}

// Cold two-stream growth rate for equal beams ±v_d, mode k:
// ω⁴ − 2 ωp² ω² − 4 k² vd² ωp² + ωp⁴ = 0  → take imag part of unstable root.
// Simplified estimate: γ ≈ (ωp / √2) · (k vd / ωp)^{2/3} for resonant cold beams (order-of-magnitude).
inline TheoryCurve twoStreamGrowthEstimate(const FormParams& params) {
	TheoryCurve curve;
	curve.name = "Two-stream estimate";
	curve.isDamping = false;
	const double k = waveNumber(params);
	const double vd = std::abs(params.driftVelocity);
	const double wp = params.plasmaFrequency;
	if (params.numSpecies < 2 || k <= 0.0 || vd <= 1e-8 || wp <= 0.0) {
		curve.note = "Need two beams with finite drift for two-stream estimate.";
		return curve;
	}
	// Exact cold equal-beam dispersion: solve for Im(ω) of
	// 1 = ωp²/2 / (ω − k vd)² + ωp²/2 / (ω + k vd)²
	// Peak growth near k vd ≈ ωp / √2. Use closed-form branch:
	const double alpha = k * vd / wp;
	if (alpha <= 0.0) {
		return curve;
	}
	// γ/ωp ≈ √( √(1 + 8 α²)/2  − (1 + 4 α²)/2 )   for unstable branch when α small-moderate
	const double inside = 0.5 * std::sqrt(1.0 + 8.0 * alpha * alpha) - 0.5 * (1.0 + 4.0 * alpha * alpha);
	if (inside <= 0.0) {
		curve.note = "Mode outside cold two-stream unstable band (estimate).";
		return curve;
	}
	curve.rate = wp * std::sqrt(inside);
	curve.valid = std::isfinite(curve.rate) && curve.rate > 0.0;
	curve.note = "Cold equal-beam two-stream growth estimate.";
	return curve;
}

inline TheoryCurve selectTheory(const FormParams& params) {
	if (params.dimension != 1) {
		TheoryCurve empty;
		empty.note = "Theory overlay currently for 1D runs.";
		return empty;
	}
	if (params.numSpecies >= 2 && std::abs(params.driftVelocity) > 1e-6 && params.thermalVelocity < 0.05) {
		return twoStreamGrowthEstimate(params);
	}
	if (params.thermalVelocity > 0.05 && std::abs(params.driftVelocity) < 0.05) {
		return landauDampingEstimate(params);
	}
	TheoryCurve empty;
	empty.note = "No matching analytic overlay for this setup.";
	return empty;
}

// Linear fit of ln|A| vs t over the middle growth/damping window.
inline GrowthFit fitExponentialRate(
	const std::vector<double>& times,
	const std::vector<double>& amplitudes,
	bool preferGrowth) {
	GrowthFit fit;
	if (times.size() < 8 || times.size() != amplitudes.size()) {
		return fit;
	}

	const size_t n = times.size();
	const size_t i0 = n / 5;
	const size_t i1 = (n * 3) / 5;
	std::vector<double> xs;
	std::vector<double> ys;
	xs.reserve(i1 - i0);
	ys.reserve(i1 - i0);
	for (size_t i = i0; i < i1; ++i) {
		if (amplitudes[i] > 1e-14) {
			xs.push_back(times[i]);
			ys.push_back(std::log(amplitudes[i]));
		}
	}
	if (xs.size() < 4) {
		return fit;
	}

	double sumX = 0.0, sumY = 0.0, sumXX = 0.0, sumXY = 0.0;
	for (size_t i = 0; i < xs.size(); ++i) {
		sumX += xs[i];
		sumY += ys[i];
		sumXX += xs[i] * xs[i];
		sumXY += xs[i] * ys[i];
	}
	const double denom = static_cast<double>(xs.size()) * sumXX - sumX * sumX;
	if (std::abs(denom) < 1e-18) {
		return fit;
	}
	const double slope = (static_cast<double>(xs.size()) * sumXY - sumX * sumY) / denom;
	const double intercept = (sumY - slope * sumX) / static_cast<double>(xs.size());

	fit.rate = slope;
	fit.amplitude0 = std::exp(intercept);
	fit.tStart = xs.front();
	fit.tEnd = xs.back();
	fit.valid = std::isfinite(slope) && std::isfinite(fit.amplitude0);
	if (preferGrowth && slope > 0.0) {
		fit.label = "Measured growth γ";
	} else if (!preferGrowth && slope < 0.0) {
		fit.label = "Measured damping γ";
	} else {
		fit.label = "Measured rate γ";
	}
	return fit;
}

} // namespace PhysicsAnalysis
