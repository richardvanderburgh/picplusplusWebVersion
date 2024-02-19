#include <chrono>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>

#include <SDL2/SDL.h>

#include <Logger.h>
#include <DataStructs.h>
#include <PICPlusPlus.h>

struct App {
    SDL_Renderer* renderer;
    SDL_Window* window;

    bool isRunning{};
} app;

void initSDL()
{
    LOG_INFO("Initializing SDL");

    int rendererFlags, windowFlags;

    rendererFlags = SDL_RENDERER_ACCELERATED;

    windowFlags = 0;

    if (SDL_Init(SDL_INIT_VIDEO) < 0)
    {
        LOG_CRITICAL("Couldn't initialize SDL: {}\n", SDL_GetError());
        return;
    }

    const int32_t SCREEN_WIDTH = 800;
    const int32_t SCREEN_HEIGHT = 600;
    app.window = SDL_CreateWindow("PIC++", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, windowFlags);

    if (!app.window)
    {
        LOG_CRITICAL("Failed to open {} x {} window: {}\n", SCREEN_WIDTH, SCREEN_HEIGHT, SDL_GetError());
        return;
    }

    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "linear");

    app.renderer = SDL_CreateRenderer(app.window, -1, rendererFlags);

    if (!app.renderer)
    {
        LOG_CRITICAL("Failed to create renderer: {}\n", SDL_GetError());
        return;
    }

    app.isRunning = true;
    LOG_INFO("SDL Initialized");
}

void doInput(void)
{
    SDL_Event event;

    while (SDL_PollEvent(&event))
    {
        switch (event.type)
        {
        case SDL_QUIT:
            app.isRunning = false;
            LOG_WARN("QUITTING");
            break;

        default:
            break;
        }
    }
}

void prepareScene(void)
{
    SDL_SetRenderDrawColor(app.renderer, 96, 128, 255, 255);
    SDL_RenderClear(app.renderer);
}

void presentScene(void)
{
    SDL_RenderPresent(app.renderer);
}

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
    Logger::Logger::Init();

    /*initSDL();

    atexit(cleanup);

    while (app.isRunning)
    {
        prepareScene();

        doInput();

        presentScene();

        SDL_Delay(16);
    }

    return 0;*/

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

	PIC_PLUS_PLUS::PICPlusPlus picPlusPlus(inputVariables);
	auto jsonResult = picPlusPlus.initialize();

	auto finish = std::chrono::high_resolution_clock::now();

	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

	std::cout << "numParticles: " << inputVariables.allSpeciesData[0].numParticles << std::endl;
	std::cout << "numTimeSteps: " << inputVariables.simulationParams.numTimeSteps << std::endl;
	std::cout << "numGrid: "	  << inputVariables.simulationParams.numGrid << std::endl;
	std::cout << "numSpecies: "   << inputVariables.simulationParams.numSpecies << std::endl;

	std::cout << "PIC++ took " << microseconds.count() << " micro secs\n";

	return 0;
}
