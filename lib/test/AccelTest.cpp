#include <gtest/gtest.h>

#include <Accel.hpp>

TEST(AccelTest, accel) {
	
	const int numSpecies = 5;
	const double gridWidth = 1.0;
	const double timeWidth = 1.0;
	const int timeStep = 0;
	const int numGrid = 1;
	const std::vector<double> particleMasses(numSpecies, 1.0);
	const std::vector<int> numSpeciesParticles(numSpecies, 1.0);
	const std::vector<std::vector<double>> particlePositions(numSpecies, std::vector<double>(1.0));
	
	double ael = 1.0;
	std::vector<double> inOutAcceleration(numGrid+1, 1.0);
	std::vector<double> inOutParticleCharges(numSpecies, 1.0);
	std::vector <std::vector<double>> inOutVelocities(numSpecies, std::vector<double>(1.0));

	accel(inOutAcceleration,
		inOutParticleCharges,
		numSpecies, gridWidth, timeWidth, timeStep,
		particleMasses,
		ael,
		numGrid,
		numSpeciesParticles,
		particlePositions,
		inOutVelocities);

	EXPECT_EQ(ael, -0.5);
	EXPECT_EQ(inOutAcceleration, std::vector<double>({-0.5, -0.5}));
	EXPECT_EQ(particlePositions, std::vector<std::vector<double>>({ { 0 }, { 0 }, { 0 }, { 0 }, { 0 } }));
	EXPECT_EQ(inOutVelocities, std::vector<std::vector<double>>({ { -0.5 }, { -0.5 }, { -0.5 }, { -0.5 }, { -0.5 } }));
}