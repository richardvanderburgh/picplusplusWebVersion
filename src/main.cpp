#include <chrono>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include <PICPlusPlus.h>

DATA_STRUCTS::InputVariables loadJSONFile(nlohmann::json config) {

	DATA_STRUCTS::InputVariables inputVariables;

	DATA_STRUCTS::SimulationParams simulationParams;
	simulationParams.numGrid = config.value("numGrid", 0);
	simulationParams.numTimeSteps = config.value("numTimeSteps", 0);
	simulationParams.spatialLength = config.value("spatialLength", 0.0);
	simulationParams.timeStepSize = config.value("timeStepSize", 0.0);
	simulationParams.numSpecies = config.value("numSpecies", 0);

	std::vector<DATA_STRUCTS::SpeciesData> allSpeciesData;

	for (const auto& speciesConfig : config["species"]) {
		DATA_STRUCTS::SpeciesData speciesData;
		speciesData.name = speciesConfig.value("name", "UnnamedSpecies");
		speciesData.numParticles = speciesConfig.value("numParticles", 0);
		speciesData.spatialPerturbationMode = speciesConfig.value("spatialPerturbationMode", 0);
		speciesData.driftVelocity = speciesConfig.value("driftVelocity", 0.0);
		speciesData.spatialPerturbationAmplitude = speciesConfig.value("spatialPerturbationAmplitude", 0.0);
		speciesData.thermalVelocity = speciesConfig.value("thermalVelocity", 0.0);
		speciesData.plasmaFrequency = speciesConfig.value("plasmaFrequency", 0);
		speciesData.chargeMassRatio = speciesConfig.value("chargeMassRatio", 0);

		std::vector<double> particlePositions(speciesData.numParticles, 0);
		std::vector<double> particleXVelocities(speciesData.numParticles, 0);

		speciesData.particlePositions = particlePositions;
		speciesData.particleXVelocities = particleXVelocities;


		allSpeciesData.push_back(speciesData);
	}
	inputVariables.simulationParams = simulationParams;
	inputVariables.allSpeciesData = allSpeciesData;

	return inputVariables;
}

int main(int argc, char* argv[]) {

	auto start = std::chrono::high_resolution_clock::now();

	if (argc != 2) {
		std::cerr << "Path: " << argv[0] << " <config.json>\n";
		return 1;
	}

	std::ifstream configFile(argv[1]);
	if (!configFile.is_open()) {
		std::cerr << "Error opening file: " << argv[1] << "\n";
		return 1;
	}

	nlohmann::json config;
	configFile >> config;

	if (config.find("species") == config.end() || !config["species"].is_array()) {
		std::cerr << "Error: 'species' array not found in the JSON file.\n";
		return 1;
	}

	DATA_STRUCTS::InputVariables inputVariables = loadJSONFile(config);
	std::cout << "numParticles: " << inputVariables.allSpeciesData[0].numParticles << std::endl;
	std::cout << "numTimeSteps: " << inputVariables.simulationParams.numTimeSteps << std::endl;
	std::cout << "numGrid: "	  << inputVariables.simulationParams.numGrid << std::endl;
	std::cout << "numSpecies: "   << inputVariables.simulationParams.numSpecies << std::endl;

	PIC_PLUS_PLUS::PICPlusPlus picPlusPlus(inputVariables);
	auto jsonResult = picPlusPlus.initialize();

	auto finish = std::chrono::high_resolution_clock::now();

	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

	std::cout << "PIC++ took " << microseconds.count() << " micro secs\n";

	return 0;
}