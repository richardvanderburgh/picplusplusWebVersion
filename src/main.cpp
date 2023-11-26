#include <chrono>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>


#include <PICPlusPlus.h>

using json = nlohmann::json;

int main(int argc, char* argv[]) {

    const double defaultL = 6.28318530717958;
    const int defaultN = 500;
    const int defaultNt = 300;
    const double defaultDt = 0.1;
    const int defaultNg = 256;
    const int defaultMode = 1;
    const double defaultV0 = 0;
    const int defaultNumSpecies = 1;
    const double defaultAmplitude = 0.02;
    const double defaultVT1 = 0.5;
    const int defaultWp1 = 1;
    const int defaultQm1 = -1;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>\n";
        return 1;
    }

    std::ifstream configFile(argv[1]);
    if (!configFile.is_open()) {
        std::cerr << "Error opening file: " << argv[1] << "\n";
        return 1;
    }

    json config;
    configFile >> config;

        // Check if the 'species' key exists
    if (config.find("species") == config.end() || !config["species"].is_array()) {
        std::cerr << "Error: 'species' array not found in the JSON file.\n";
        return 1;
    }

    // Extract species data
    const int numSpecies = config.value("numSpecies", 0);
    std::vector<DATA_STRUCTS::SpeciesData> allSpeciesData ;

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

    const double spatialLength = config.value("spatialLength", 0.0);
    const int numTimeSteps = config.value("numTimeSteps", 0);
    const double timeStepSize = config.value("timeStepSize", 0.0);
    const int numGrid = config.value("numGrid", 0);

    auto start = std::chrono::high_resolution_clock::now();

    PIC_PLUS_PLUS::PICPlusPlus init(
        spatialLength,
        numTimeSteps,
        timeStepSize,
        numGrid,
        numSpecies,
        allSpeciesData);

	auto jsonResult = init.initialize();

    auto finish = std::chrono::high_resolution_clock::now();

    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    std::cout << "numParticles: " << allSpeciesData[0].numParticles << std::endl;
    std::cout << "numTimeSteps: " << numTimeSteps << std::endl;
    std::cout << "numGrid: "      << numGrid << std::endl;
    std::cout << "numSpecies: "   << numSpecies << std::endl;

    std::cout << "PIC++ took " << microseconds.count() << " micro secs\n";

	return 0;
}
