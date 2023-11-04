#ifndef SET_RHO_HPP
#define SET_RHO_HPP

#include <vector>

void setRho(int species, int ng, double dx, 
	std::vector<int> N, 
	std::vector<double> qdx, 
	std::vector<double>& rho, 
	std::vector<std::vector<double>>& x, 
	std::vector<std::vector<double>>& rho0, 
	std::vector<std::vector<double>>& rhos, int iw) {

	//double dxi;
	const double dxi = 1.0 / dx;
	const int xn = ng;

	// If it is the first group of particles, then clear out rho.
	if (species == 0) {
		for (int j = 1; j < ng + 1; j++) {
			rho[j] = rho0[0][j];
		}
		rho[0] = 0;
	}

	// Add on fixed neutralizing charge density
	// (not needed when all species are mobile - but harmless.)
	for (int j = 0; j < rhos[0].size(); j++) {
		rhos[species][j] = 0;
	}
	for (int j = 1; j <= ng; j++) {
		rho0[species][j] = rho0[species][j] - rhos[species][j];
		rho[j] = rho[j] - rhos[species][j];
	}

	// NGP
	if (iw == 1) {
		for (int i = 1; i <= N[species]; i++) {
			x[species][i] = x[species][i] * dxi;
			if (x[species][i] < 0) {
				x[species][i] = x[species][i] + xn;
			}
			if (x[species][i] > xn) {
				x[species][i] = x[species][i] - xn;
			}
			int64_t j = static_cast<int64_t>(floor((x[species][i]) + 1 + 0.5));
			rho[j] = rho[j] + qdx[species];
		}
	}

	// Linear
	else if (iw == 2) {
		for (int i = 0; i < N[species]; i++) {
			x[species][i] = x[species][i] * dxi;
			if (x[species][i] < 0) {
				x[species][i] = x[species][i] + xn;
			}
			if (x[species][i] > xn) {
				x[species][i] = x[species][i] - xn;
			}
			int64_t j = static_cast<int64_t>(floor(x[species][i]));
			//jdata(i, species) = j;
			double drho = qdx[species] * (x[species][i] - j);
			//drhodata(i, species) = drho;
			rho[j] = rho[j] - drho + qdx[species];
			rho[j + 1] = rho[j + 1] + drho;
		}
	}
}


void move(int nsp, std::vector<double>& rho, std::vector<std::vector<double>> rho0, std::vector<double> qdx, std::vector<int> N, std::vector <std::vector<double>>& x, std::vector <std::vector<double>>& vx, int ng) {
	//MOVE - Advances position one time step and accumulates charge density.
	for (int species = 0; species < nsp; species++) {
		// Clear out old charge density.
		for (int j = 1; j <= ng; j++) {
			rho[j] = rho0[species][j];
		}
		rho[0] = 0;
	}

	for (int species = 0; species < nsp; species++) {
		for (int i = 0; i < N[species]; i++) {

			x[species][i] += vx[species][i];

			if (x[species][i] < 0)
				x[species][i] += ng;

			if (x[species][i] >= ng)
				x[species][i] -= ng;


			int64_t j = static_cast<int64_t>(floor(x[species][i]));
			double drho = qdx[species] * (x[species][i] - j);
			rho[j] = rho[j] - drho + qdx[species];
			rho[j + 1] = rho[j + 1] + drho;
		}
	}

}
#endif // !SET_RHO_HPP