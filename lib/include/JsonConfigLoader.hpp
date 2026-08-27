#pragma once

#include <nlohmann/json.hpp>

#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include <DataStructs.h>
#include "Utils.hpp"

inline DATA_STRUCTS::InputVariables loadJSONFile(const nlohmann::json& config) {
	DATA_STRUCTS::InputVariables inputVariables;

	DATA_STRUCTS::SimulationParams simulationParams;
	simulationParams.numGrid = config.value("numGrid", 0);
	simulationParams.numTimeSteps = config.value("numTimeSteps", 0);
	simulationParams.spatialLength = config.value("spatialLength", 0.0);
	simulationParams.timeStepSize = config.value("timeStepSize", 0.0);
	simulationParams.numSpecies = config.value("numSpecies", 0);
	simulationParams.framePeriod = config.value("framePeriod", 1);

	std::vector<DATA_STRUCTS::SpeciesData> allSpeciesData;

	for (const auto& speciesConfig : config["species"]) {
		DATA_STRUCTS::SpeciesData speciesData;
		speciesData.name = speciesConfig.value("name", "UnnamedSpecies");
		speciesData.numParticles = speciesConfig.value("numParticles", 0);
		speciesData.spatialPerturbationMode = speciesConfig.value("spatialPerturbationMode", 0);
		speciesData.driftVelocity = speciesConfig.value("driftVelocity", 0.0);
		speciesData.spatialPerturbationAmplitude = speciesConfig.value("spatialPerturbationAmplitude", 0.0);
		speciesData.spatialPerturbationWaveform = speciesConfig.value("spatialPerturbationWaveform", "cos");
		speciesData.thermalVelocity = speciesConfig.value("thermalVelocity", 0.0);
		speciesData.plasmaFrequency = speciesConfig.value("plasmaFrequency", 0);
		speciesData.chargeMassRatio = speciesConfig.value("chargeMassRatio", 0);

		speciesData.particlePositions = std::vector<double>(speciesData.numParticles, 0);
		speciesData.particleXVelocities = std::vector<double>(speciesData.numParticles, 0);

		allSpeciesData.push_back(speciesData);
	}

	inputVariables.simulationParams = simulationParams;
	inputVariables.allSpeciesData = allSpeciesData;

	return inputVariables;
}

inline std::optional<DATA_STRUCTS::InputVariables> loadJSONFileFromPath(const std::string& path, std::string& error) {
	std::ifstream configFile(path);
	if (!configFile.is_open()) {
		error = "Error opening file: " + path;
		return std::nullopt;
	}

	nlohmann::json config;
	try {
		configFile >> config;
	} catch (const std::exception& ex) {
		error = std::string("Invalid JSON: ") + ex.what();
		return std::nullopt;
	}

	if (config.find("species") == config.end() || !config["species"].is_array()) {
		error = "Error: 'species' array not found in the JSON file.";
		return std::nullopt;
	}

	return loadJSONFile(config);
}
