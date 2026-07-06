#pragma once 

#include <cmath>
#include <cstddef>
#include <vector>

void setRho(int species,
	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	const std::vector<double>& qdx,
	std::vector<double>& rho,
	std::vector<std::vector<double>>& rho0,
	std::vector<std::vector<double>>& rhos) {

	const double dxi = 1.0 / simulationParams.gridStepSize;
	const int xn = simulationParams.numGrid;

	// If it is the first group of particles, then clear out rho.
	if (species == 0) {
		for (int j = 1; j < simulationParams.numGrid + 1; j++) {
			rho[j] = rho0[0][j];
		}
		rho[0] = 0;
	}

	// Add on fixed neutralizing charge density
	// (not needed when all species are mobile - but harmless.)
	for (int j = 0; j < rhos[0].size(); j++) {
		rhos[species][j] = 0;
	}
	for (int j = 1; j <= simulationParams.numGrid; j++) {
		rho0[species][j] = rho0[species][j] - rhos[species][j];
		rho[j] = rho[j] - rhos[species][j];
	}

	for (int i = 0; i < allSpeciesData[species].numParticles; i++) {
		allSpeciesData[species].particlePositions[i] = allSpeciesData[species].particlePositions[i] * dxi;
		if (allSpeciesData[species].particlePositions[i] < 0) {
			allSpeciesData[species].particlePositions[i] = allSpeciesData[species].particlePositions[i] + xn;
		}
		if (allSpeciesData[species].particlePositions[i] > xn) {
			allSpeciesData[species].particlePositions[i] = allSpeciesData[species].particlePositions[i] - xn;
		}
		int64_t j = static_cast<int64_t>(floor(allSpeciesData[species].particlePositions[i]));
		double drho = qdx[species] * (allSpeciesData[species].particlePositions[i] - j);
		rho[j] = rho[j] - drho + qdx[species];
		rho[j + 1] = rho[j + 1] + drho;
	}
}

void move(
	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	std::vector<double>& rho, 
	const std::vector<std::vector<double>>& rho0, 
	const std::vector<double>& qdx) {

	for (int species = 0; species < simulationParams.numSpecies; species++) {
		// Clear out old charge density.
		for (int j = 1; j <= simulationParams.numGrid; j++) {
			rho[j] = rho0[species][j];
		}
		rho[0] = 0;
	}

	for (int species = 0; species < simulationParams.numSpecies; species++) {
		std::vector<double>& positions = allSpeciesData[species].particlePositions;
		const std::vector<double>& velocities = allSpeciesData[species].particleXVelocities;
		const int numParticles = allSpeciesData[species].numParticles;
		const double q = qdx[species];
		const double gridLength = static_cast<double>(simulationParams.numGrid);

		// The position update is independent per particle, but the charge
		// deposition is a scatter (rho[j], rho[j+1]) that races when threads
		// touch the same cell. Each thread accumulates into a private density
		// buffer and merges once under a critical section, so the result is
		// race-free and does not serialize the inner loop with atomics.
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			std::vector<double> localRho(rho.size(), 0.0);

#ifdef _OPENMP
#pragma omp for nowait schedule(static)
#endif
			for (int i = 0; i < numParticles; i++) {
				double position = positions[i] + velocities[i];

				if (position < 0.0)
					position += gridLength;

				if (position >= gridLength)
					position -= gridLength;

				positions[i] = position;

				const double gridPosition = std::floor(position);
				const size_t j = static_cast<size_t>(gridPosition);
				const double drho = q * (position - gridPosition);
				localRho[j] += q - drho;
				localRho[j + 1] += drho;
			}

#ifdef _OPENMP
#pragma omp critical
#endif
			{
				for (size_t k = 0; k < rho.size(); ++k) {
					rho[k] += localRho[k];
				}
			}
		}
	}

}