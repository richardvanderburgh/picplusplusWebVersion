#include <chrono>
#include <fstream>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <PICPlusPlus.h>
#include "JsonConfigLoader.hpp"

int main(int argc, char* argv[]) {

	auto start = std::chrono::high_resolution_clock::now();

	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <config.json> [output.json]\n";
		return 1;
	}

	std::string error;
	const auto inputVariables = loadJSONFileFromPath(argv[1], error);
	if (!inputVariables.has_value()) {
		std::cerr << error << "\n";
		return 1;
	}

	if (const auto validationError = validateSimulationParams(inputVariables->simulationParams)) {
		std::cerr << "Invalid simulation parameters: " << *validationError << "\n";
		return 1;
	}

	std::cout << "numParticles: " << inputVariables->allSpeciesData[0].numParticles << std::endl;
	std::cout << "numTimeSteps: " << inputVariables->simulationParams.numTimeSteps << std::endl;
	std::cout << "numGrid: "	  << inputVariables->simulationParams.numGrid << std::endl;
	std::cout << "numSpecies: "   << inputVariables->simulationParams.numSpecies << std::endl;
	std::cout << "dimension: "    << inputVariables->simulationParams.dimension << std::endl;
	if (inputVariables->simulationParams.dimension == 3) {
		const auto& p = inputVariables->simulationParams;
		const int ny = p.numGridY > 0 ? p.numGridY : p.numGrid;
		const int nz = p.numGridZ > 0 ? p.numGridZ : p.numGrid;
		std::cout << "numGridY: " << ny << std::endl;
		std::cout << "numGridZ: " << nz << std::endl;
	}
#ifdef _OPENMP
	std::cout << "OpenMP threads: " << omp_get_max_threads() << std::endl;
#else
	std::cout << "OpenMP threads: 1 (serial build)" << std::endl;
#endif

	PIC_PLUS_PLUS::PICPlusPlus picPlusPlus(*inputVariables);
	auto jsonResult = picPlusPlus.initialize();

	if (!jsonResult.has_value()) {
		std::cerr << "Simulation failed to produce results.\n";
		return 1;
	}

	if (argc >= 3) {
		std::ofstream outputFile(argv[2]);
		if (!outputFile.is_open()) {
			std::cerr << "Error opening output file: " << argv[2] << "\n";
			return 1;
		}
		outputFile << jsonResult->dump(2);
	}

	auto finish = std::chrono::high_resolution_clock::now();

	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

	std::cout << "PIC++ took " << microseconds.count() << " micro secs\n";

	return 0;
}
