#ifndef ACCEL_HPP
#define ACCEL_HPP

#include <cmath>
#include <vector>

inline void accel(

	DATA_STRUCTS::SimulationParams simulationParams,
	
	std::vector<double>& inOutAcceleration, 
	std::vector<double> particleCharges,
	int t,
	const std::vector<double>& particleMasses, 
	double& ael,
	const std::vector<int>& N, 
	const std::vector <std::vector<double>>& particlePositions,
	std::vector <std::vector<double>>& inOutVelocities) {

	const double dxdt = simulationParams.gridStepSize / simulationParams.timeStepSize;
	for (int species = 0; species < simulationParams.numSpecies; species++) {

		// This looks weird but it's because the first time step is -1/2 a step

		if (t == 0)
			particleCharges[species] = -0.5 * particleCharges[species];

		const double ae = (particleCharges[species] / particleMasses[species]) * (simulationParams.timeStepSize / dxdt);

		//  renormalizes acceleration if need be.
		if (ae != ael) {
			const double tem = ae / ael;
			for (int j = 0; j <= simulationParams.numGrid; j++) {
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