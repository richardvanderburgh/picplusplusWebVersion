#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>

#include <nlohmann/json.hpp>

#include <Fields.hpp>
#include <PICPlusPlus.h>

namespace {

	DATA_STRUCTS::InputVariables buildTwoStreamInput()
	{
		DATA_STRUCTS::SimulationParams simulationParams;
		simulationParams.numGrid = 32;
		simulationParams.numTimeSteps = 500;
		simulationParams.spatialLength = 6.28318530717958;
		simulationParams.timeStepSize = 0.2;
		simulationParams.numSpecies = 2;

		auto makeSpecies = [](double driftVelocity) {
			DATA_STRUCTS::SpeciesData species;
			species.numParticles = 500;
			species.spatialPerturbationMode = 1;
			species.driftVelocity = driftVelocity;
			species.spatialPerturbationAmplitude = 0.001;
			species.thermalVelocity = 0;
			species.plasmaFrequency = 1;
			species.chargeMassRatio = -1;
			species.particlePositions = std::vector<double>(species.numParticles, 0.0);
			species.particleXVelocities = std::vector<double>(species.numParticles, 0.0);
			return species;
		};

		DATA_STRUCTS::InputVariables inputVariables;
		inputVariables.simulationParams = simulationParams;
		inputVariables.allSpeciesData = { makeSpecies(1.0), makeSpecies(-1.0) };
		return inputVariables;
	}

	double sumKineticEnergy(const nlohmann::json& result)
	{
		double total = 0.0;
		for (const auto& speciesKe : result["ke"]) {
			total += speciesKe.back().get<double>();
		}
		return total;
	}

} // namespace

TEST(ValidationTest, TwoStreamInstabilityDevelops)
{
	const auto inputVariables = buildTwoStreamInput();
	PIC_PLUS_PLUS::PICPlusPlus simulation(inputVariables);

	const auto result = simulation.initialize();
	ASSERT_TRUE(result.has_value());

	const auto& ese = result->at("ese");
	const double initialFieldEnergy = ese.front().get<double>();
	const double finalFieldEnergy = ese.back().get<double>();

	EXPECT_GT(finalFieldEnergy, initialFieldEnergy * 5.0)
		<< "Two-stream instability should amplify the initial field perturbation";
}

TEST(ValidationTest, TwoStreamEnergyBudget)
{
	const auto inputVariables = buildTwoStreamInput();
	PIC_PLUS_PLUS::PICPlusPlus simulation(inputVariables);

	const auto result = simulation.initialize();
	ASSERT_TRUE(result.has_value());

	const double initialTotal = sumKineticEnergy(*result) + result->at("ese").front().get<double>();
	const double finalTotal = sumKineticEnergy(*result) + result->at("ese").back().get<double>();
	const double relativeDrift = std::abs(finalTotal - initialTotal) / initialTotal;

	EXPECT_LT(relativeDrift, 0.15)
		<< "Total energy should remain within 15% over the two-stream run";
}

TEST(ValidationTest, FieldsSolverUsesConfiguredGridSize)
{
	DATA_STRUCTS::SimulationParams simulationParams;
	simulationParams.numGrid = 256;
	simulationParams.spatialLength = 6.28318530717958;
	simulationParams.gridStepSize = simulationParams.spatialLength / simulationParams.numGrid;

	std::vector<double> chargeDensity(simulationParams.numGrid + 1, 0.0);
	chargeDensity[1] = 1.0;
	chargeDensity[simulationParams.numGrid - 1] = -1.0;

	std::vector<std::vector<double>> electricField(1, std::vector<double>(simulationParams.numGrid + 1, 0.0));
	std::vector<double> acceleration(simulationParams.numGrid + 1, 0.0);
	double ael = 1.0;

	fields(simulationParams, chargeDensity, electricField, 0, acceleration, ael);

	double maxField = 0.0;
	for (double value : electricField[0]) {
		maxField = std::max(maxField, std::abs(value));
	}

	EXPECT_GT(maxField, 1e-6)
		<< "Field solver should produce a non-zero electric field at numGrid=256";
}
