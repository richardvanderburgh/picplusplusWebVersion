#include <chrono>
#include <iostream>

#include <SDL2/SDL.h>

#include <Logger.h>

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

int main(int argc, char* argv[]) {
    Logger::Logger::Init();

    initSDL();

    //atexit(cleanup);

    while (app.isRunning)
    {
        prepareScene();

        doInput();

        presentScene();

        SDL_Delay(16);
    }

    return 0;


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

    // Initialize variables with default values or command line arguments
    const double spatialLength = argc > 1 ? std::stod(argv[1]) : defaultL;
    const int numParticles = argc > 2 ? atoi(argv[2]) : defaultN;
    const int numTimeSteps = argc > 3 ? atoi(argv[3]) : defaultNt;
    const double timeStepSize = argc > 4 ? std::stod(argv[4]) : defaultDt;
    const int numGrid = argc > 5 ? atoi(argv[5]) : defaultNg;
    const int spatialPerturbationMode = argc > 6 ? std::stoi(argv[6]) : defaultMode;
    const double driftVelocity = argc > 7 ? std::stod(argv[7]) : defaultV0;
    const int numSpecies = argc > 8 ? std::stoi(argv[8]) : defaultNumSpecies;
    const double spatialPerturbationAmplitude = argc > 9 ? std::stod(argv[9]) : defaultAmplitude;
    const double thermalVelocity = argc > 10 ? std::stod(argv[10]) : defaultVT1;
    const int plasmaFrequency = argc > 11 ? std::stoi(argv[11]) : defaultWp1;
    const int chargeMassRatio = argc > 12 ? std::stoi(argv[12]) : defaultQm1;

    // 6.28318530717958 5 3 0.1 32 1 1 2 0.001 0.0 1.0 -1.0

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
    std::cout << "Finished in " << microseconds.count() << " micro secs\n";

	return 0;
}
