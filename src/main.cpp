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

    // Check if a JSON file is provided as a command line argument
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>\n";
        return 1;
    }

    // Read the JSON file
    std::ifstream configFile(argv[1]);
    if (!configFile.is_open()) {
        std::cerr << "Error opening file: " << argv[1] << "\n";
        return 1;
    }

    json config;
    configFile >> config;


    // Initialize variables with default values or values from the JSON file
    const double spatialLength = config.value("spatialLength", defaultL);
    const int numParticles = config.value("numParticles", defaultN);
    const int numTimeSteps = config.value("numTimeSteps", defaultNt);
    const double timeStepSize = config.value("timeStepSize", defaultDt);
    const int numGrid = config.value("numGrid", defaultNg);
    const int spatialPerturbationMode = config.value("spatialPerturbationMode", defaultMode);
    const double driftVelocity = config.value("driftVelocity", defaultV0);
    const int numSpecies = config.value("numSpecies", defaultNumSpecies);
    const double spatialPerturbationAmplitude = config.value("spatialPerturbationAmplitude", defaultAmplitude);
    const double thermalVelocity = config.value("thermalVelocity", defaultVT1);
    const int plasmaFrequency = config.value("plasmaFrequency", defaultWp1);
    const int chargeMassRatio = config.value("chargeMassRatio", defaultQm1);

    auto start = std::chrono::high_resolution_clock::now();

	PIC_PLUS_PLUS::PICPlusPlus init(spatialLength,
        numParticles,
        numTimeSteps,
        timeStepSize,
        numGrid,
        spatialPerturbationMode,
        driftVelocity,
        numSpecies,
        spatialPerturbationAmplitude,
        thermalVelocity,
        plasmaFrequency,
        chargeMassRatio);

	auto jsonResult = init.initialize();

    auto finish = std::chrono::high_resolution_clock::now();

    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    std::cout << "numParticles: " << numParticles << std::endl;
    std::cout << "numTimeSteps: " << numTimeSteps << std::endl;
    std::cout << "numGrid: "      << numGrid << std::endl;
    std::cout << "numSpecies: "   << numSpecies << std::endl;

    std::cout << "PIC++ took " << microseconds.count() << " micro secs\n";

	return 0;
}
