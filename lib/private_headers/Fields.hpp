#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "complex.hpp"
#include "fft.hpp"
#include <DataStructs.h>

inline void fields(
	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<double>& inOutChargeDensity,
	std::vector<std::vector<double>>& inOutElectricField,
	const int timeStep,
	std::vector<double>& inOutAcceleration,
	double& ael) {

	const int numGrid = simulationParams.numGrid;
	std::vector<double> electricPotential(numGrid + 1, 0.0);

	inOutChargeDensity[0] += inOutChargeDensity[numGrid];
	inOutChargeDensity[numGrid] = inOutChargeDensity[0];

	std::vector<complex> complexChargeDensity(numGrid + 1);
	std::vector<complex> complexChargeDensityK(numGrid);
	std::vector<complex> complexPotential(numGrid + 1);
	std::vector<complex> complexPotentialK(numGrid + 1);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numGrid >= 1024)
#endif
	for (int i = 0; i < numGrid; i++) {
		complexChargeDensity[i] = inOutChargeDensity[i];
	}

	if (!CFFT::Forward(complexChargeDensity.data(), complexChargeDensityK.data(),
			static_cast<unsigned int>(numGrid))) {
		throw std::runtime_error("FFT forward failed: numGrid must be a power of two");
	}

	const double lengthScale = simulationParams.spatialLength;
	const double twoPiOverL = 2.0 * std::numbers::pi / lengthScale;
	const int halfGrid = numGrid / 2;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numGrid >= 1024)
#endif
	for (int k = 0; k < numGrid + 1; k++) {
		int mode = 0;
		if (k == 0) {
			mode = 1;
		} else if (k <= halfGrid) {
			mode = k;
		} else if (k > halfGrid) {
			mode = k - numGrid;
		}

		if (mode == 0) {
			complexPotentialK[k] = complex(0.0, 0.0);
			continue;
		}

		const double kSquared = std::pow(twoPiOverL * static_cast<double>(mode), 2.0);
		complexPotentialK[k] = complexChargeDensityK[k] / (-kSquared);
	}

	if (!CFFT::Inverse(complexPotentialK.data(), complexPotential.data(),
			static_cast<unsigned int>(numGrid))) {
		throw std::runtime_error("FFT inverse failed: numGrid must be a power of two");
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numGrid >= 1024)
#endif
	for (int i = 0; i < numGrid + 1; i++) {
		electricPotential[i] = -complexPotential[i].re();
	}

	electricPotential[numGrid] = electricPotential[0];

	const double invTwoDx = 1.0 / (2.0 * simulationParams.gridStepSize);
	std::vector<double>& electricField = inOutElectricField[timeStep];

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(numGrid >= 1024)
#endif
	for (int j = 1; j < numGrid; j++) {
		electricField[j] = (electricPotential[j - 1] - electricPotential[j + 1]) * invTwoDx;
	}

	electricField[0] = (electricPotential[numGrid - 1] - electricPotential[1]) * invTwoDx;
	electricField[numGrid] = electricField[0];

	ael = 1;
	inOutAcceleration = electricField;
}
#endif
