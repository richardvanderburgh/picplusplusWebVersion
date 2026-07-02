#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <cmath>
#include <numbers>
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

	for (int i = 0; i < numGrid; i++) {
		complexChargeDensity[i] = inOutChargeDensity[i];
	}
	CFFT::Forward(complexChargeDensity.data(), complexChargeDensityK.data(), numGrid);

	for (int k = 0; k < numGrid + 1; k++) {
		int ii = 0;
		if (k == 0) {
			ii = 1;
		}
		else if (k <= numGrid / 2 && k >= 0) {
			ii = k;
		}
		else if (k > numGrid / 2) {
			ii = k - numGrid;
		}

		if (ii == 0)
			break;

		complexPotentialK[k] = complexChargeDensityK[k] / -std::pow(2.0 * std::numbers::pi * ii / simulationParams.spatialLength, 2.0);
	}

	CFFT::Inverse(complexPotentialK.data(), complexPotential.data(), numGrid);

	for (int i = 0; i < numGrid + 1; i++) {
		electricPotential[i] = -complexPotential[i].re();
	}

	electricPotential[numGrid] = electricPotential[0];

	for (int j = 1; j < numGrid; j++) {
		inOutElectricField[timeStep][j] = (electricPotential[j - 1] - electricPotential[j + 1]) / (2.0 * simulationParams.gridStepSize);
	}

	inOutElectricField[timeStep][0] = (electricPotential[numGrid - 1] - electricPotential[1]) / (2.0 * simulationParams.gridStepSize);
	inOutElectricField[timeStep][numGrid] = inOutElectricField[timeStep][0];

	ael = 1;
	inOutAcceleration = inOutElectricField[timeStep];
}
#endif
