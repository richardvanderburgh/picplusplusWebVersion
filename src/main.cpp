#include <chrono>
#include <iostream>

#include <init.h>

int main(int argc, char* argv[]) {

	PIC_PLUS_PLUS::Init init;

#if 0
    double L = atoi(argv[1]);
	int N = atoi(argv[2]);
	int nt = atoi(argv[3]);
	double dt = std::stod(argv[4]);
    int ng = atoi(argv[5]);
    int mode = std::stod(argv[6]);
    double V0 = std::stod(argv[7]);
    int numSpecies = std::stod(argv[8]);
    double amplitude = std::stod(argv[9]);
    double VT1 = std::stod(argv[10]);
    int wp1 = std::stod(argv[11]);
    int qm1 = std::stod(argv[12]);
#else
    const double L = 6.28318530717958;
    const int N = 500;
    const int nt = 300;
    const double dt = 0.1;
    const int ng = 256;
    const int mode = 1;
    const double V0 = 0;
    const int numSpecies = 1;
    const double amplitude = 0.02;
    const double VT1 = 0.5;
    const int wp1 = 1;
    const int qm1 = -1;
#endif

    // 6.28318530717958 5 3 0.1 32 1 1 2 0.001 0.0 1.0 -1.0

    auto start = std::chrono::high_resolution_clock::now();
	auto jsonResult = init.initialize(L, N, nt, dt, ng, mode, V0, numSpecies, amplitude, VT1, wp1, qm1 );
    auto finish = std::chrono::high_resolution_clock::now();

    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
    std::cout << "Finished in " << microseconds.count() << " micro secs\n";

	return 0;
}