#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <cmath>
#include <numbers>

#include "complex.hpp"
#include "fft.hpp"


inline void fields(
	DATA_STRUCTS::SimulationParams simulationParams,
	std::vector<double>& inOutChargeDensity,
	std::vector<std::vector<double>>& inOutElectricField, 
	const int timeStep,  
	std::vector<double>& inOutAcceleration, 
	double& ael) {
	
	const int numGrid = 32;
	std::vector<double> electricPotential(numGrid + 1, 0.0);
	double l = simulationParams.spatialLength;

	inOutChargeDensity[0] += inOutChargeDensity[numGrid];
	inOutChargeDensity[numGrid] = inOutChargeDensity[0];

	complex complexChargeDensity[numGrid + 1];
	complex complexChargeDensityK[numGrid];

	complex complexPotential[numGrid + 1];
	complex complexPotentialK[numGrid + 1];

	//fftw_complex complexRho[ng + 1];
	//fftw_complex complexRhok[ng];   // Output array
	//fftw_complex complexPhik[ng + 1];
	//fftw_complex complexPhi[ng + 1];
	//fftw_plan forwardPlan;

	for (int i = 0; i < numGrid; i++) { complexChargeDensity[i] = inOutChargeDensity[i]; }
	CFFT::Forward(complexChargeDensity, complexChargeDensityK, numGrid);

	//for (int i = 0; i < ng + 1; i++) {
		//complexRho[i][0] = rho[i];
		//complexRho[i][1] = 0;
		//complexPhik[i][0] = 0;
		//complexPhik[i][1] = 0;
		//complexPhi[i][0] = 0;
		//complexPhi[i][1] = 0;
	//}

	//forwardPlan = fftw_plan_dft_1d(ng, complexRho, complexRhok, FFTW_FORWARD,  FFTW_ESTIMATE);
	//fftw_execute(forwardPlan);
	//fftw_destroy_plan(forwardPlan);

	for (int k = 0; k < numGrid + 1; k++) {
		int ii = 0;
		if (k == 0) {
			ii = 1;
		}
		else if (k <= numGrid / 2 && k >= 0) {
			ii = k;
		}
		else if (k > numGrid / 2) {
			ii = k - numGrid;
		}

		if (ii == 0)
			break;

		complexPotentialK[k] = complexChargeDensityK[k] / -std::pow(2.0 * std::numbers::pi * ii / simulationParams.spatialLength, 2.0);

		//complexPhik[k][0] = complexRhok[k][0] / -std::pow(2.0 * M_PI * ii / L, 2.0);
		//complexPhik[k][1] = complexRhok[k][1] / -std::pow(2.0 * M_PI * ii / L, 2.0);

	}

	//complex complexPhi[ng];
	CFFT::Inverse(complexPotentialK, complexPotential, numGrid);

	//fftw_plan inversePlan;
	//inversePlan = fftw_plan_dft_1d(ng, complexPhik, complexPhi, FFTW_BACKWARD, FFTW_ESTIMATE);
	//fftw_execute(inversePlan);

	//fftw_destroy_plan(inversePlan);
	//fftw_free(complexRhok);

	//for (int i = 0; i < ng + 1; i++) {
	//	double real = complexPhi[i][0];
	//	double imag = complexPhi[i][1];

	//	phi[i] = -real;
	//}

	for (int i = 0; i < numGrid + 1; i++) {
		electricPotential[i] = -complexPotential[i].re();
	}

	//fftw_free(complexPhi);


	electricPotential[numGrid] = electricPotential[0];

		for (int j = 1; j < numGrid; j++) {
			inOutElectricField[timeStep][j] = (electricPotential[j - 1] - electricPotential[j + 1]) / (2.0 * simulationParams.gridStepSize);
		}

		inOutElectricField[timeStep][0] = (electricPotential[numGrid - 1] - electricPotential[1]) / (2.0 * simulationParams.gridStepSize);
		inOutElectricField[timeStep][numGrid] = inOutElectricField[timeStep][0];


	ael = 1;
	for (int i = 0; i <= numGrid; i++) {
		inOutAcceleration[i] = inOutElectricField[timeStep][i];
	}
}
#endif