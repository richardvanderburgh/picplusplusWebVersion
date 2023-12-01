#include <gtest/gtest.h>

#include <Fields.hpp>
#include "DataStructs.h"

TEST(FieldsTest, fields) {
	std::vector<double> chargeDensity(33, 1.0);
	double spatialLength = 1.0;
	double gridSpacing = 1.0;
	std::vector<std::vector<double>> electricField(2, std::vector<double>(33, 1.0));
	int timeStep = 1;
	const int numGrid = 1;
	std::vector<double> gridAccelleration(33, 0.0);
	double ael = 1.0;

	DATA_STRUCTS::SimulationParams simulationParams;
	simulationParams.spatialLength = spatialLength;
	simulationParams.numGrid = numGrid;
	fields(simulationParams, chargeDensity, electricField, timeStep, gridAccelleration, ael);

	//fields(rho, L, dx, E, t, constNg, a, ael);

	//EXPECT_EQ(rho, std::vector<double>());
	//EXPECT_EQ(E, std::vector<std::vector<double>>());
	//EXPECT_EQ(a, std::vector<double>());
	EXPECT_EQ(ael, 1.0);
}
