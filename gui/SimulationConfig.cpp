#include "SimulationConfig.h"

#include <fstream>

#include "Utils.hpp"

namespace {

std::string demoManifestPath(const std::string& repoRoot) {
	return repoRoot + "/inputFiles/demo/manifest.json";
}

std::string demoConfigPath(const std::string& repoRoot, const std::string& file) {
	return repoRoot + "/inputFiles/demo/" + file;
}

} // namespace

std::vector<DemoEntry> SimulationConfig::loadDemoManifest(const std::string& repoRoot, std::string& error) {
	std::ifstream manifestFile(demoManifestPath(repoRoot));
	if (!manifestFile.is_open()) {
		error = "Could not open demo manifest.";
		return {};
	}

	nlohmann::json manifest;
	try {
		manifestFile >> manifest;
	} catch (const std::exception& ex) {
		error = std::string("Invalid demo manifest: ") + ex.what();
		return {};
	}

	std::vector<DemoEntry> demos;
	for (const auto& entry : manifest) {
		DemoEntry demo;
		demo.id = entry.value("id", "");
		demo.file = entry.value("file", "");
		demo.title = entry.value("title", demo.id);
		demo.category = entry.value("category", "");
		demo.description = entry.value("description", "");
		if (entry.contains("lookFor") && entry["lookFor"].is_array()) {
			for (const auto& item : entry["lookFor"]) {
				demo.lookFor.push_back(item.get<std::string>());
			}
		}
		if (!demo.id.empty()) {
			demos.push_back(demo);
		}
	}
	return demos;
}

FormParams SimulationConfig::formParamsFromJsonConfig(const nlohmann::json& config) {
	FormParams params;
	if (!config.contains("species") || config["species"].empty()) {
		return params;
	}

	const auto& species0 = config["species"][0];
	params.spatialLength = config.value("spatialLength", params.spatialLength);
	params.numParticles = species0.value("numParticles", params.numParticles);
	params.timeSteps = config.value("numTimeSteps", params.timeSteps);
	params.timeStepSize = config.value("timeStepSize", params.timeStepSize);
	params.numGrid = config.value("numGrid", params.numGrid);
	params.spatialPerturbationMode = species0.value("spatialPerturbationMode", params.spatialPerturbationMode);
	params.driftVelocity = species0.value("driftVelocity", params.driftVelocity);
	params.numSpecies = config.value("numSpecies", params.numSpecies);
	params.spatialPerturbationAmplitude = species0.value("spatialPerturbationAmplitude", params.spatialPerturbationAmplitude);
	params.thermalVelocity = species0.value("thermalVelocity", params.thermalVelocity);
	params.plasmaFrequency = species0.value("plasmaFrequency", params.plasmaFrequency);
	params.chargeMassRatio = species0.value("chargeMassRatio", params.chargeMassRatio);
	params.spatialPerturbationWaveform = species0.value("spatialPerturbationWaveform", params.spatialPerturbationWaveform);
	params.framePeriod = config.value("framePeriod", params.framePeriod);
	if (params.framePeriod < 1) {
		params.framePeriod = 5;
	}
	return params;
}

std::optional<FormParams> SimulationConfig::formParamsFromDemo(
	const std::string& repoRoot, const std::string& demoId, std::string& error) {
	const auto demos = loadDemoManifest(repoRoot, error);
	if (!error.empty()) {
		return std::nullopt;
	}

	const DemoEntry* match = nullptr;
	for (const auto& demo : demos) {
		if (demo.id == demoId) {
			match = &demo;
			break;
		}
	}
	if (match == nullptr) {
		error = "Unknown demo: " + demoId;
		return std::nullopt;
	}

	std::ifstream demoFile(demoConfigPath(repoRoot, match->file));
	if (!demoFile.is_open()) {
		error = "Could not open demo config: " + match->file;
		return std::nullopt;
	}

	nlohmann::json config;
	try {
		demoFile >> config;
	} catch (const std::exception& ex) {
		error = std::string("Invalid demo JSON: ") + ex.what();
		return std::nullopt;
	}

	return formParamsFromJsonConfig(config);
}

nlohmann::json SimulationConfig::buildConfig(const FormParams& params) {
	nlohmann::json speciesTemplate = {
		{"numParticles", params.numParticles},
		{"spatialPerturbationMode", params.spatialPerturbationMode},
		{"spatialPerturbationAmplitude", params.spatialPerturbationAmplitude},
		{"spatialPerturbationWaveform", params.spatialPerturbationWaveform},
		{"thermalVelocity", params.thermalVelocity},
		{"plasmaFrequency", params.plasmaFrequency},
		{"chargeMassRatio", params.chargeMassRatio},
	};

	nlohmann::json species = nlohmann::json::array();
	if (params.numSpecies == 1) {
		auto s = speciesTemplate;
		s["name"] = "Species1";
		s["driftVelocity"] = params.driftVelocity;
		species.push_back(s);
	} else {
		auto s1 = speciesTemplate;
		s1["name"] = "Species1";
		s1["driftVelocity"] = params.driftVelocity;
		species.push_back(s1);

		auto s2 = speciesTemplate;
		s2["name"] = "Species2";
		s2["driftVelocity"] = -params.driftVelocity;
		species.push_back(s2);
	}

	const int framePeriod = params.framePeriod < 1 ? 1 : params.framePeriod;

	return {
		{"species", species},
		{"spatialLength", params.spatialLength},
		{"numTimeSteps", params.timeSteps},
		{"timeStepSize", params.timeStepSize},
		{"numGrid", params.numGrid},
		{"numSpecies", static_cast<int>(species.size())},
		{"framePeriod", framePeriod},
	};
}

std::optional<std::string> SimulationConfig::validateFormParams(const FormParams& params) {
	if (!isPowerOfTwo(params.numGrid)) {
		return "Grid cells must be a power of 2 (e.g. 32, 64, 128, 256).";
	}
	if (params.numParticles < 1) {
		return "Need at least 1 particle per species.";
	}
	if (params.timeSteps < 1) {
		return "Need at least 1 time step.";
	}
	if (params.timeStepSize <= 0.0) {
		return "Time step Δt must be positive.";
	}
	if (params.numSpecies < 1 || params.numSpecies > 2) {
		return "Species count must be 1 or 2.";
	}
	return std::nullopt;
}
