#ifndef UTILS_HPP
#define UTILS_HPP

#include <optional>
#include <string>
#include <vector>

#include "DataStructs.h"

inline bool isPowerOfTwo(int value) {
	return value > 0 && (value & (value - 1)) == 0;
}

// Returns an error message when simulation parameters are invalid.
inline std::optional<std::string> validateSimulationParams(
	const DATA_STRUCTS::SimulationParams& params) {
	if (params.numGrid < 1) {
		return "numGrid must be positive (got " + std::to_string(params.numGrid) + ")";
	}
	if (!isPowerOfTwo(params.numGrid)) {
		return "numGrid must be a power of two for the FFT field solver (got "
			+ std::to_string(params.numGrid) + ")";
	}
	if (params.numTimeSteps < 0) {
		return "numTimeSteps must be non-negative";
	}
	if (params.numSpecies < 1) {
		return "numSpecies must be at least 1";
	}
	if (params.spatialLength <= 0.0) {
		return "spatialLength must be positive";
	}
	if (params.timeStepSize <= 0.0) {
		return "timeStepSize must be positive";
	}
	return std::nullopt;
}

inline std::vector<double> linspace(double start, double end, int num_points) {
	std::vector<double> result(num_points);
	const double step = (end - start) / (num_points - 1);
	for (int i = 0; i < num_points; ++i) {
		result[i] = start + i * step;
	}
	return result;
}
#endif // !UTILS_HPP
