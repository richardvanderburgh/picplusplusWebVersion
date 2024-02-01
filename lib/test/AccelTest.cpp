#include <gtest/gtest.h>

#include <Accel.hpp>
#include "DataStructs.h"

TEST(AccelTest, accel) {
	DATA_STRUCTS::SimulationParams simulationParams;

	simulationParams.numSpecies = 5;

	const std::vector<double> particleMasses(simulationParams.numSpecies, 1.0);
	const std::vector<double> particleCharges(simulationParams.numSpecies, 1.0);
	const std::vector<int> numSpeciesParticles(simulationParams.numSpecies, 1.0);

	const std::vector<std::vector<double>> particlePositions(simulationParams.numSpecies, std::vector<double>(1.0));
	std::vector <std::vector<double>> particleXVelocities(simulationParams.numSpecies, std::vector<double>(1.0));
	
	simulationParams.gridStepSize = 1;
	simulationParams.timeStepSize = 1;
	simulationParams.numGrid = 1;

	std::vector<DATA_STRUCTS::SpeciesData> allSpeciesData(simulationParams.numSpecies);
	for (int i = 0; i < simulationParams.numSpecies; i++) {
		allSpeciesData[i].particleMass = particleMasses[i];
		allSpeciesData[i].particleCharge = particleCharges[i];

		allSpeciesData[i].numParticles = 1;
		allSpeciesData[i].particlePositions = particlePositions[i];
		allSpeciesData[i].particleXVelocities = particleXVelocities[i];
	}

	std::vector<double> inOutAcceleration(simulationParams.numGrid + 1, 1.0);
	const int timeStep = 0;
	double ael = 1.0;

	accel(simulationParams, allSpeciesData, inOutAcceleration, timeStep, ael);

	EXPECT_EQ(ael, -0.5);
	EXPECT_EQ(inOutAcceleration, std::vector<double>({-0.5, -0.5}));

	for (int i = 0; i < simulationParams.numSpecies; i++) {
		EXPECT_EQ(allSpeciesData[i].particlePositions[0], 0);
		EXPECT_EQ(allSpeciesData[i].particleXVelocities[0], -0.5);
	}
}