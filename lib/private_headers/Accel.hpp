#ifndef ACCEL_HPP
#define ACCEL_HPP

#include <cmath>
#include <vector>

#include "Boris.hpp"
#include "DataStructs.h"

inline void accel(

	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	std::vector<double>& inOutAcceleration, 
	int t,
	double& ael) {

	const double dxdt = simulationParams.gridStepSize / simulationParams.timeStepSize;
	const bool magnetized = simulationParams.hasMagneticField();
	const double bx = simulationParams.magneticFieldX;
	const double by = simulationParams.magneticFieldY;
	const double bz = simulationParams.magneticFieldZ;
	const double dt = simulationParams.timeStepSize;
	const double pushDt = (t == 0) ? -0.5 * dt : dt;

	for (int species = 0; species < simulationParams.numSpecies; species++) {

		const double qm = allSpeciesData[species].particleCharge / allSpeciesData[species].particleMass;
		double ae = qm * (pushDt / dxdt);

		// Renormalize the shared acceleration grid when species q/m changes.
		// Grid stores ae·E after rescale so the unmagnetized gather can add it directly.
		if (ae != ael) {
			const double tem = ae / ael;
			for (int j = 0; j <= simulationParams.numGrid; j++) {
				inOutAcceleration[j] *= tem;
			}
			ael = ae;
		}

		const auto& acceleration = inOutAcceleration;
		const std::vector<double>& positions = allSpeciesData[species].particlePositions;
		std::vector<double>& velocities = allSpeciesData[species].particleXVelocities;
		const int numParticles = allSpeciesData[species].numParticles;

		if (!magnetized) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (int i = 0; i < numParticles; ++i) {
				const double gridPosition = std::floor(positions[i]);
				const size_t j = static_cast<size_t>(gridPosition);
				velocities[i] = velocities[i] + acceleration[j] +
					(positions[i] - gridPosition) * (acceleration[j + 1] - acceleration[j]);
			}
			continue;
		}

		// 1D3V Boris: Ex from the 1D mesh; Ey = Ez = 0. Transverse velocities
		// use the same dx/dt scale as vx (standard 1D3V convention).
		std::vector<double>& vy = allSpeciesData[species].particleYVelocities;
		std::vector<double>& vz = allSpeciesData[species].particleZVelocities;
		if (vy.size() != static_cast<size_t>(numParticles) || vz.size() != static_cast<size_t>(numParticles)) {
			continue;
		}

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (int i = 0; i < numParticles; ++i) {
			const double gridPosition = std::floor(positions[i]);
			const size_t j = static_cast<size_t>(gridPosition);
			// acceleration[] already holds ae·Ex after species rescale.
			const double aeEx = acceleration[j] +
				(positions[i] - gridPosition) * (acceleration[j + 1] - acceleration[j]);
			const double exSample = (ae != 0.0) ? (aeEx / ae) : 0.0;
			borisPushCodeUnits(
				velocities[i], vy[i], vz[i],
				exSample, 0.0, 0.0,
				ae, ae, ae,
				qm, pushDt,
				dxdt, dxdt, dxdt,
				bx, by, bz);
		}
	}
}
#endif
