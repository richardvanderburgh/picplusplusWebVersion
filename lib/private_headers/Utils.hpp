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
	if (params.dimension != 1 && params.dimension != 3) {
		return "dimension must be 1 or 3 (got " + std::to_string(params.dimension) + ")";
	}
	if (params.dimension == 3) {
		const int ny = params.numGridY > 0 ? params.numGridY : params.numGrid;
		const int nz = params.numGridZ > 0 ? params.numGridZ : params.numGrid;
		if (!isPowerOfTwo(ny)) {
			return "numGridY must be a power of two for the 3D FFT field solver (got "
				+ std::to_string(ny) + ")";
		}
		if (!isPowerOfTwo(nz)) {
			return "numGridZ must be a power of two for the 3D FFT field solver (got "
				+ std::to_string(nz) + ")";
		}
		const double ly = params.spatialLengthY > 0.0 ? params.spatialLengthY : params.spatialLength;
		const double lz = params.spatialLengthZ > 0.0 ? params.spatialLengthZ : params.spatialLength;
		if (ly <= 0.0) {
			return "spatialLengthY must be positive";
		}
		if (lz <= 0.0) {
			return "spatialLengthZ must be positive";
		}
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
