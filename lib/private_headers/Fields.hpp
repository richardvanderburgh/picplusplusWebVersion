#ifndef FIELDS_HPP
#define FIELDS_HPP

#include <cmath>
#include <numbers>

#include "complex.hpp"
#include "fft.hpp"


inline void fields(std::vector<double>& rho,
	double L, double dx, 
	std::vector<std::vector<double>>& E, 
	int t, const int constNg, 
	std::vector<double>& a, 
	double& ael) {
	
	const int ng = 32;
	std::vector<double> phi(ng + 1, 0.0);
	double l = L;

	rho[0] += rho[ng];
	rho[ng] = rho[0];

	complex complexRho[ng + 1];
	complex complexRhok[ng];

	complex complexPhik[ng + 1];
	complex complexPhi[ng + 1];

	//fftw_complex complexRho[ng + 1];
	//fftw_complex complexRhok[ng];   // Output array
	//fftw_complex complexPhik[ng + 1];
	//fftw_complex complexPhi[ng + 1];
	//fftw_plan forwardPlan;

	for (int i = 0; i < ng; i++) { complexRho[i] = rho[i]; }
	CFFT::Forward(complexRho, complexRhok, ng);

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

	for (int k = 0; k < ng + 1; k++) {
		int ii = 0;
		if (k == 0) {
			ii = 1;
		}
		else if (k <= ng / 2 && k >= 0) {
			ii = k;
		}
		else if (k > ng / 2) {
			ii = k - ng;
		}

		if (ii == 0)
			break;

		complexPhik[k] = complexRhok[k] / -std::pow(2.0 * std::numbers::pi * ii / L, 2.0);

		//complexPhik[k][0] = complexRhok[k][0] / -std::pow(2.0 * M_PI * ii / L, 2.0);
		//complexPhik[k][1] = complexRhok[k][1] / -std::pow(2.0 * M_PI * ii / L, 2.0);

	}

	//complex complexPhi[ng];
	CFFT::Inverse(complexPhik, complexPhi, ng);

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

	for (int i = 0; i < ng + 1; i++) {
		phi[i] = -complexPhi[i].re();
	}

	//fftw_free(complexPhi);


	phi[ng] = phi[0];

		for (int j = 1; j < ng; j++) {
			E[t][j] = (phi[j - 1] - phi[j + 1]) / (2.0 * dx);
		}

		E[t][0] = (phi[ng - 1] - phi[1]) / (2.0 * dx);
		E[t][ng] = E[t][0];


	ael = 1;
	for (int i = 0; i <= ng; i++) {
		a[i] = E[t][i];
	}
}
#endif