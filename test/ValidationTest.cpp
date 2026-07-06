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

	DATA_STRUCTS::InputVariables buildLandauDampingInput()
	{
		DATA_STRUCTS::SimulationParams simulationParams;
		simulationParams.numGrid = 256;
		simulationParams.numTimeSteps = 500;
		simulationParams.spatialLength = 6.28318530717958;
		simulationParams.timeStepSize = 0.1;
		simulationParams.numSpecies = 1;
		simulationParams.framePeriod = 0;

		DATA_STRUCTS::SpeciesData species;
		species.numParticles = 8000;
		species.spatialPerturbationMode = 1;
		species.spatialPerturbationWaveform = "sin";
		species.driftVelocity = 0;
		species.spatialPerturbationAmplitude = 0.02;
		species.thermalVelocity = 0.4;
		species.plasmaFrequency = 1;
		species.chargeMassRatio = -1;
		species.particlePositions = std::vector<double>(species.numParticles, 0.0);
		species.particleXVelocities = std::vector<double>(species.numParticles, 0.0);

		DATA_STRUCTS::InputVariables inputVariables;
		inputVariables.simulationParams = simulationParams;
		inputVariables.allSpeciesData = { species };
		return inputVariables;
	}

	DATA_STRUCTS::InputVariables buildColdPlasmaWaveInput()
	{
		auto inputVariables = buildLandauDampingInput();
		inputVariables.allSpeciesData[0].thermalVelocity = 0;
		return inputVariables;
	}

	double peakFieldEnergy(const nlohmann::json& result)
	{
		const auto& ese = result.at("ese");
		double peak = 0.0;
		for (const auto& value : ese) {
			peak = std::max(peak, value.get<double>());
		}
		return peak;
	}

	double sumKineticEnergy(const nlohmann::json& result)
	{
		double total = 0.0;
		for (const auto& speciesKe : result["ke"]) {
			total += speciesKe.back().get<double>();
		}
		return total;
	}

	double meanFieldEnergy(const nlohmann::json& result, int beginStep, int endStep)
	{
		const auto& ese = result.at("ese");
		double sum = 0.0;
		int count = 0;
		for (int step = beginStep; step < endStep && step < static_cast<int>(ese.size()); ++step) {
			sum += ese[step].get<double>();
			++count;
		}
		return count > 0 ? sum / static_cast<double>(count) : 0.0;
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

TEST(ValidationTest, LandauDampingReducesFieldEnergy)
{
	const auto inputVariables = buildLandauDampingInput();
	PIC_PLUS_PLUS::PICPlusPlus simulation(inputVariables);

	const auto result = simulation.initialize();
	ASSERT_TRUE(result.has_value());

	const auto& ese = result->at("ese");
	const double initialFieldEnergy = ese.front().get<double>();
	const double finalFieldEnergy = ese.back().get<double>();
	const double peakFieldEnergyValue = peakFieldEnergy(*result);

	EXPECT_GT(peakFieldEnergyValue, initialFieldEnergy * 1.05)
		<< "Landau wave should briefly exchange energy into the field";
	EXPECT_LT(finalFieldEnergy, peakFieldEnergyValue * 0.55)
		<< "Landau damping should reduce field energy well below its peak";
	EXPECT_LT(finalFieldEnergy, initialFieldEnergy * 0.85)
		<< "Net field energy should decay below the initial perturbation level";
}

TEST(ValidationTest, LandauEnergyBudget)
{
	const auto inputVariables = buildLandauDampingInput();
	PIC_PLUS_PLUS::PICPlusPlus simulation(inputVariables);

	const auto result = simulation.initialize();
	ASSERT_TRUE(result.has_value());

	const double initialTotal = sumKineticEnergy(*result) + result->at("ese").front().get<double>();
	const double finalTotal = sumKineticEnergy(*result) + result->at("ese").back().get<double>();
	const double relativeDrift = std::abs(finalTotal - initialTotal) / initialTotal;

	EXPECT_LT(relativeDrift, 0.15)
		<< "Total energy should remain within 15% over the Landau damping run";
}

TEST(ValidationTest, WarmPlasmaDampsMoreThanCold)
{
	auto warmInput = buildLandauDampingInput();
	auto coldInput = buildColdPlasmaWaveInput();

	const auto warmResult = PIC_PLUS_PLUS::PICPlusPlus(warmInput).initialize();
	const auto coldResult = PIC_PLUS_PLUS::PICPlusPlus(coldInput).initialize();
	ASSERT_TRUE(warmResult.has_value());
	ASSERT_TRUE(coldResult.has_value());

	const double warmLateToEarly = meanFieldEnergy(*warmResult, 400, 500) / meanFieldEnergy(*warmResult, 50, 150);
	const double coldLateToEarly = meanFieldEnergy(*coldResult, 400, 500) / meanFieldEnergy(*coldResult, 50, 150);

	EXPECT_LT(warmLateToEarly, 0.85)
		<< "Warm plasma should show a clear decay in mean field energy over the run";
	EXPECT_GT(coldLateToEarly, warmLateToEarly + 0.1)
		<< "Cold plasma should damp less than the warm Landau-damped case";
}

TEST(ValidationTest, ColdPlasmaWaveRemainsUndamped)
{
	const auto inputVariables = buildColdPlasmaWaveInput();
	const auto result = PIC_PLUS_PLUS::PICPlusPlus(inputVariables).initialize();
	ASSERT_TRUE(result.has_value());

	const double lateToEarly = meanFieldEnergy(*result, 400, 500) / meanFieldEnergy(*result, 50, 150);
	EXPECT_GT(lateToEarly, 0.90)
		<< "Cold plasma should not show sustained Landau damping of the launched wave";
	EXPECT_LT(lateToEarly, 1.20)
		<< "Field energy should stay near its initial level, not grow without bound";
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

TEST(ValidationTest, RejectsNonPowerOfTwoGrid)
{
	auto inputVariables = buildTwoStreamInput();
	inputVariables.simulationParams.numGrid = 30;

	PIC_PLUS_PLUS::PICPlusPlus simulation(inputVariables);
	const auto result = simulation.initialize();

	EXPECT_FALSE(result.has_value())
		<< "numGrid=30 should be rejected because the FFT solver requires a power of two";
}

TEST(ValidationTest, ThermalVelocityIsReproducible)
{
	auto inputA = buildTwoStreamInput();
	auto inputB = buildTwoStreamInput();
	inputA.allSpeciesData[0].thermalVelocity = 0.1;
	inputA.allSpeciesData[1].thermalVelocity = 0.1;
	inputB.allSpeciesData[0].thermalVelocity = 0.1;
	inputB.allSpeciesData[1].thermalVelocity = 0.1;
	inputA.simulationParams.numTimeSteps = 10;
	inputB.simulationParams.numTimeSteps = 10;

	const auto resultA = PIC_PLUS_PLUS::PICPlusPlus(inputA).initialize();
	const auto resultB = PIC_PLUS_PLUS::PICPlusPlus(inputB).initialize();

	ASSERT_TRUE(resultA.has_value());
	ASSERT_TRUE(resultB.has_value());
	EXPECT_DOUBLE_EQ(sumKineticEnergy(*resultA), sumKineticEnergy(*resultB));
}
