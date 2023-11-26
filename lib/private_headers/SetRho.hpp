#pragma once 

#include <vector>

void setRho(int species,
	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	const std::vector<double>& qdx,
	std::vector<double>& rho,
	std::vector<std::vector<double>>& rho0,
	std::vector<std::vector<double>>& rhos) {

	//double dxi;
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
		//jdata(i, species) = j;
		double drho = qdx[species] * (allSpeciesData[species].particlePositions[i] - j);
		//drhodata(i, species) = drho;
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
		for (int i = 0; i < allSpeciesData[species].numParticles; i++) {

			allSpeciesData[species].particlePositions[i] += allSpeciesData[species].particleXVelocities[i];

			if (allSpeciesData[species].particlePositions[i] < 0)
				allSpeciesData[species].particlePositions[i] += simulationParams.numGrid;

			if (allSpeciesData[species].particlePositions[i] >= simulationParams.numGrid)
				allSpeciesData[species].particlePositions[i] -= simulationParams.numGrid;


			int64_t j = static_cast<int64_t>(floor(allSpeciesData[species].particlePositions[i]));
			double drho = qdx[species] * (allSpeciesData[species].particlePositions[i] - j);
			rho[j] = rho[j] - drho + qdx[species];
			rho[j + 1] = rho[j + 1] + drho;
		}
	}

}