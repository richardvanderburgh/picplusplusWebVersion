#ifndef ACCEL_HPP
#define ACCEL_HPP

#include <cmath>
#include <vector>

inline void accel(
	
	std::vector<double>& inOutAcceleration, 
	std::vector<double> particleCharges, 
	int nsp, double dx, double dt, int t,
	const std::vector<double>& particleMasses, 
	double& ael, int numGrid, 
	const std::vector<int>& N, 
	const std::vector <std::vector<double>>& particlePositions,
	std::vector <std::vector<double>>& inOutVelocities) {

	const double dxdt = dx / dt;
	for (int species = 0; species < nsp; species++) {

		// This looks weird but it's because the first time step is -1/2 a step

		if (t == 0)
			particleCharges[species] = -0.5 * particleCharges[species];

		const double ae = (particleCharges[species] / particleMasses[species]) * (dt / dxdt);

		//  renormalizes acceleration if need be.
		if (ae != ael) {
			const double tem = ae / ael;
			for (int j = 0; j <= numGrid; j++) {
				inOutAcceleration[j] *= tem;
			}
			ael = ae;
		}

		for (int i = 0; i < N[species]; ++i) {
			const int64_t j = static_cast<int64_t>(floor(particlePositions[species][i]));
			inOutVelocities[species][i] = inOutVelocities[species][i] + inOutAcceleration[j] + 
				(particlePositions[species][i] - j) * (inOutAcceleration[j + 1] - inOutAcceleration[j]);
		}
	}
}
#endif