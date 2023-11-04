
#include <init.h>

int main(int argc, char* argv[]) {


	PIC_PLUS_PLUS::Init init;
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

    // 6.28318530717958 5 3 0.1 32 1 1 2 0.001 0.0 1.0 -1.0

	bool success = init.initialize(L, N, nt, dt, ng, mode, V0, numSpecies, amplitude, VT1, wp1, qm1 );

	return 0;
}