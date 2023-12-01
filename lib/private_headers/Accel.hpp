#ifndef ACCEL_HPP
#define ACCEL_HPP

#include <cmath>
#include <vector>

#include "DataStructs.h"

inline void accel(

	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	std::vector<double>& inOutAcceleration, 
	int t,
	double& ael) {

	const double dxdt = simulationParams.gridStepSize / simulationParams.timeStepSize;
	for (int species = 0; species < simulationParams.numSpecies; species++) {

		// This looks weird but it's because the first time step is -1/2 a step
		double ae = (allSpeciesData[species].particleCharge / allSpeciesData[species].particleMass) * (simulationParams.timeStepSize / dxdt);

		if (t == 0)
			ae = -0.5 * ae;


		//  renormalizes acceleration if need be.
		if (ae != ael) {
			const double tem = ae / ael;
			for (int j = 0; j <= simulationParams.numGrid; j++) {
				inOutAcceleration[j] *= tem;
			}
			ael = ae;
		}

		for (int i = 0; i < allSpeciesData[species].numParticles; ++i) {
			const int64_t j = static_cast<int64_t>(floor(allSpeciesData[species].particlePositions[i]));
			allSpeciesData[species].particleXVelocities[i] = allSpeciesData[species].particleXVelocities[i] + inOutAcceleration[j] +
				(allSpeciesData[species].particlePositions[i] - j) * (inOutAcceleration[j + 1] - inOutAcceleration[j]);
		}
	}
}
#endif