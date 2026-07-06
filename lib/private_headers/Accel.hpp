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

		// Gather-only particle push: every particle reads the shared acceleration
		// grid and writes its own velocity, so this loop is embarrassingly parallel
		// (no races, no reduction).
		const auto& acceleration = inOutAcceleration;
		const std::vector<double>& positions = allSpeciesData[species].particlePositions;
		std::vector<double>& velocities = allSpeciesData[species].particleXVelocities;
		const int numParticles = allSpeciesData[species].numParticles;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (int i = 0; i < numParticles; ++i) {
			const double gridPosition = std::floor(positions[i]);
			const size_t j = static_cast<size_t>(gridPosition);
			velocities[i] = velocities[i] + acceleration[j] +
				(positions[i] - gridPosition) * (acceleration[j + 1] - acceleration[j]);
		}
	}
}
#endif