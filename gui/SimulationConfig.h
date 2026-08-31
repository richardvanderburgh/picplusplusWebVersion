#pragma once

#include <QMetaType>

#include <nlohmann/json.hpp>

#include <optional>
#include <string>
#include <vector>

struct FormParams {
	int dimension = 1;
	double spatialLength = 6.28318530717958;
	double spatialLengthY = 6.28318530717958;
	double spatialLengthZ = 6.28318530717958;
	int numParticles = 500;
	int timeSteps = 500;
	double timeStepSize = 0.2;
	int numGrid = 32;
	int numGridY = 8;
	int numGridZ = 8;
	int spatialPerturbationMode = 1;
	double driftVelocity = 1.0;
	int numSpecies = 2;
	double spatialPerturbationAmplitude = 0.001;
	double thermalVelocity = 0.0;
	double plasmaFrequency = 1.0;
	double chargeMassRatio = -1.0;
	std::string spatialPerturbationWaveform = "cos";
	int framePeriod = 5;
};

Q_DECLARE_METATYPE(FormParams)

struct DemoEntry {
	std::string id;
	std::string file;
	std::string title;
	std::string category;
	std::string description;
	std::vector<std::string> lookFor;
};

struct LessonStep {
	std::string demoId;
	std::string prompt;
};

struct LessonEntry {
	std::string id;
	std::string title;
	std::string description;
	std::vector<LessonStep> steps;
};

class SimulationConfig {
public:
	static std::vector<DemoEntry> loadDemoManifest(const std::string& repoRoot, std::string& error);
	static std::vector<LessonEntry> loadLessons(const std::string& repoRoot, std::string& error);
	static std::optional<FormParams> formParamsFromDemo(const std::string& repoRoot, const std::string& demoId, std::string& error);
	static std::optional<nlohmann::json> loadDemoJson(const std::string& repoRoot, const std::string& demoId, std::string& error);
	static nlohmann::json buildConfig(const FormParams& params);
	static std::optional<std::string> validateFormParams(const FormParams& params);
	static FormParams formParamsFromJsonConfig(const nlohmann::json& config);
};
