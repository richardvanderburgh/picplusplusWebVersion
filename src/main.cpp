#include <chrono>
#include <iostream>

#include <init.h>

int main(int argc, char* argv[]) {

	PIC_PLUS_PLUS::PICPlusPlus init;

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
    const double L = argc > 1 ? std::stod(argv[1]) : defaultL;
    const int N = argc > 2 ? atoi(argv[2]) : defaultN;
    const int nt = argc > 3 ? atoi(argv[3]) : defaultNt;
    const double dt = argc > 4 ? std::stod(argv[4]) : defaultDt;
    const int ng = argc > 5 ? atoi(argv[5]) : defaultNg;
    const int mode = argc > 6 ? std::stoi(argv[6]) : defaultMode;
    const double V0 = argc > 7 ? std::stod(argv[7]) : defaultV0;
    const int numSpecies = argc > 8 ? std::stoi(argv[8]) : defaultNumSpecies;
    const double amplitude = argc > 9 ? std::stod(argv[9]) : defaultAmplitude;
    const double VT1 = argc > 10 ? std::stod(argv[10]) : defaultVT1;
    const int wp1 = argc > 11 ? std::stoi(argv[11]) : defaultWp1;
    const int qm1 = argc > 12 ? std::stoi(argv[12]) : defaultQm1;

    // 6.28318530717958 5 3 0.1 32 1 1 2 0.001 0.0 1.0 -1.0

    auto start = std::chrono::high_resolution_clock::now();
	auto jsonResult = init.initialize(L, N, nt, dt, ng, mode, V0, numSpecies, amplitude, VT1, wp1, qm1 );
    auto finish = std::chrono::high_resolution_clock::now();

    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    std::cout << "Finished in " << microseconds.count() << " micro secs\n";

	return 0;
}