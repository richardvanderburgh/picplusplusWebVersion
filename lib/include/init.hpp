#pragma once
#define _USE_MATH_DEFINES

#include <fft.hpp>
#include <fftw3.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <math.h>
// #include <corecrt_math_defines.h>

#include <iostream>
#include <thread>
#include <vector>

#include "Accel.hpp"
#include "Fields.hpp"
#include "SetRho.hpp"
#include "Utils.hpp"

class Init {
public:

	struct Particle {
		double position;
		double velocity;
		int species;
		int id;
	};

	struct Frame {
		std::vector<Particle> particles;
		std::vector<double> electricField;
		int frameNumber = 0;
	};

	struct PicData {
		std::vector<Frame> frames;
	};

	PicData mPicData;

	bool initialize(int N1, int nt, double dt, int MODE, double V0, int nsp, double amplitude, double VT1) {
		
		// Set initial input values
		
		std::string example = "Electron - Electron Stream";
		// Input Variables
		double L = 6.28318530717958; // Physical length of system in meters
		//int nsp = 3; // Number of particle species
		//int nt = 5;//600; // Number of time steps
		//double dt = .1; // Time step in seconds
		int epsi = 1; // 1 over epsilon naught(F / m) Epsilon normalized to 1
		const int ng = 32; // Number of spatial grid points - only change the power of 2
		int iw = 2; // 1 for nearest grid point, 2 for linear
		int a1 = 0; // Smoothing factor
		int a2 = 0;
		int e0 = 0; // Add uniform electric field e0* cos(w0 * time)
		int w0 = 0;

		// Plotting Intervals
		int irho = 0; // nt / 10 + 4;
		int iphi = 0; // nt / 10 + 4;
		int iE = 0; // nt / 10 + 4;
		int mplot = 3;
		int ixvx = 1; // 60
		int ivxvy = 0;
		int ifvx = 0;// nt / 10 + 4;
		int nplot = 30; // ? ?

		int wp1 = 1;
		int wc1 = 0;
		int qm1 = -1;
		//int VT1	 = 0;
		int VT2	 = 0;
		int NV2	 = 0;
		int NLG = 1;
		//int v01	 = ;
		int pch1 = 0;

		// //Species Input Variables
		std::vector<int> N(nsp);          // Number of simulation particles
		std::vector<double> wp(nsp);         // Plasma Frequency
		std::vector<double> wc(nsp);         // Cyclotron frequency
		std::vector<double> qm(nsp);         // q / m charge to mass ratio(C / kg)
		std::vector<double> vt1(nsp);        // RMS thermal velocity for random velocities
		std::vector<double> vt2(nsp);        // RMS thermal velocity for ordered velocities
		std::vector<int> nv2(nsp);	      
		std::vector<int> nlg(nsp) ;       // Number of loading groups
		std::vector<double> v0 = { 

			V0, -V0, 0, 2*V0, -2*V0};	  // Drift velocity

		std::vector<int> pch(nsp);		  // species pitch angle
		int distribution = 1;		      // Distribution 0 = Cold 1 = Two - Stream

		// Perturbation
		//int MODE = 1;
		//double amplitude = 0.001;
		int vPerturb = 0;
		int THETAX = 0;
		int THETAV = 0;

		//////////////////////////////////////////////////////////
		// INIT - Initial values for each species
		std::vector<int>    mode(nsp);
		std::vector<double> x1(nsp);
		std::vector<int>    v1(nsp);
		std::vector<int>    thetav(nsp);
		std::vector<int>    thetax(nsp);

		int npt = 0;
		for (int i = 0; i < nsp; i++) {
			N[i] = N1;
			npt += N[i];

			wp[i] = wp1;
			wc[i] = wc1;
			qm[i] = qm1;
			vt1[i] = VT1;
			vt2[i] = VT2;
			nv2[i] = NV2;
			nlg[i] = NLG;
			mode[i] = MODE;
			x1[i] = amplitude;
			v1[i] = vPerturb;
			thetax[i] = THETAX;
			thetav[i] = THETAV;
		}

		double step_size = L / ng;
		std::vector<double> gridx;

		for (double x = 0.0; x <= L; x += step_size) {
			gridx.push_back(x);
		}

		// --Field Initialization --

		std::vector<std::vector<double>> E(nt + 1, std::vector<double>(ng + 1, 0.0));

		for (int i = 0; i < ng + 1; i++) {
			for (int j = 0; j < nt + 1; j++) {
				E[j][i] = 0.0;
			}
		}

		int maxN = *std::max_element(N.begin(), N.end());

		// Energies
		// Total time for plot
		std::vector<std::vector<double>> esem(nt + 1, std::vector<double>(mplot + 1, 0.0));
		std::vector<std::vector<double>> ke(nsp, std::vector<double>(nt +1, 0.0));
		std::vector<std::vector<double>> p(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> de(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> therme(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> v_old(maxN, std::vector<double>(nsp));
		std::vector<std::vector<double>> vy_old(maxN, std::vector<double>(nsp));
		std::vector<double> ESE(nt + 1, 0.0);
		std::vector<double> esestot(nt + 1, 0.0);
		std::vector<double> te(nt + 1, 0.0);

		// positions and velocities
		std::vector<std::vector<double>> x(nsp, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> vx(nsp, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> vy(nsp, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> vz(nsp, std::vector<double>(maxN, 0.0));

		std::vector<double> ddx(nsp);
		std::vector<double> q(nsp, 0.0);
		std::vector<double> m(nsp);
		std::vector<double> T(nsp);

		for (int species = 0; species < nsp; species++)
		{
			q[species] = L * wp[species] * wp[species] / (epsi * N[species] * qm[species]);
			m[species] = q[species] / qm[species];
			T[species] = tan(-wc[species] * dt / 2.);
			ddx[species] = L / N[species];

			double ngr = N[species] / nlg[species];
			double nm = N[species] * m[species];
			double lg = L / nlg[species];

			// Set evenly spaced charge distribution distribution
			// remember that ddx is the width of the charge cloud
			for (int I = 1; I < N[species] + 1; I++) {
				x[species][I - 1] = (I - 0.5) * ddx[species];
			}

			for (int I = 0; I < N[species]; I++) {
				vx[species][I] = v0[species];
			}

			// Load ordered velocities in vx ("quiet start", or at least subdued).
			// Is set up for Maxwellian*v**nv2, but can do any smooth distribution.
			// Hereafter, ngr is preferably a power of 2.
			// First store indefinite integral of distribution function in x array.
			// Use midpoint rule -simple and quite accurate.

			int i1 = 0;

			if (vt2[species] != 0) {
				int vmax = 5 * vt2[species];
				double dv = 2 * vmax / (N[species] - 1);
				double vvnv2 = 1;
				x[species][1] = 0;

				for (int ith = 2; ith <= N[species]; ith++) {
					double vv = ((ith - 1.5) * dv - vmax) / vt2[species];

					if (nv2[species] != 0) {
						vvnv2 = pow(vv, nv2[species]);
					}

					double fv = vvnv2 * exp(-0.5 * pow(vv, 2));
					int i1 = ith - 1 + 1;
					x[species][i1] = x[species][i1 - 1] + std::max(fv, 0.0);
				}

				// For evenly spaced(half - integer multiples) values of the integral,
				// find corresponding velocity by inverse linear interpolation.

				double df = x[species][i1] / ngr;
				int i1 = 1;
				int j = 1;

				for (int ith = 0; ith < ngr; ith++) {
					double fv = (ith + 1 - 0.5) * df;

					while (fv >= x[species][j + 1]) {
						j++;
					}

					double vv = dv * (j - 1 - 1 + (fv - x[species][j - 1]) / (x[species][j] - x[species][j - 1])) - vmax;
					vx[species][i1] = vx[species][i1] + vv;
					i1++;
				}

				// For ordered velocities, scramble positions to reduce correlations.
				// Scrambling is done by bit - reversed counter - compare sorter in cpft.
				// xs = .000, .100, .010, .110, .001, .101, .011, .111, .0001.. (binary fractions)

				double xs = 0;

				for (int ith = 1; ith <= ngr; ith++) {
					int i1 = ith - 1 + 1;
					x[i1][species] = xs * lg + 0.5 * ddx[species];
					double xsi = 1.0;

					while (xs >= 0) {
						xsi = 0.5 * xsi;
						xs = xs - xsi;
					}

					xs = xs + 2.0 * xsi;
					i1 = ngr + 1 - 1;
				}
			}

			if (wc[species] != 0) {
				for (int ith = 1; ith <= ngr; ++ith) {
					int i1 = ith - 1 + 1;
					double vv = vx[species][i1];
					double theta = 2 * M_PI * x[species][i1] / lg;
					vx[species][i1] = vv * cos(theta);
					vy[species][i1] = vv * sin(theta);
				}
			}

			if (nlg[species] != 1) {
				int j = ngr + 1;
				double xs = 0;
				for (int ith = j; ith <= N[species]; ++ith) {
					xs = xs + lg;
					for (int jj = 1; jj <= ngr; ++jj) {
						int i1 = jj - 1 + 1;
						int i2 = i1 + ith - 1;
						x[species][i2] = x[species][i1] + xs;
						vx[species][i2] = vx[species][i1];
						if (wc[species] != 0) {
							vy[species][i2] = vy[species][i1];
						}
					}
				}
			}

			// Add random Maxwellian.
			if (vt1[species] != 0) {
				for (int I = 0; I < N[species]; ++I) {
					double rm = 0;
					for (int ith = 0; ith < 12; ++ith) {
						rm += (double)rand() / RAND_MAX;
					}
					rm -= 6;

					if (wc[species] != 0) {
						vy[species][I] += vt1[species] * rm;
					}
					vx[species][I] += vt1[species] * rm;
				}
			}

			for (int a = 0; a < N[species]; ++a) {
				double theta = 2 * M_PI * mode[species] * x[species][a] / L;
				x[species][a] = x[species][a] + x1[species] * cos(theta + thetax[species]);
				vx[species][a] = vx[species][a] + v1[species] * sin(theta + thetav[species]);

				if (x[species][a] >= L) {
					x[species][a] = x[species][a] - L;
				}
				if (x[species][a] < 0) {
					x[species][a] = x[species][a] + L;
				}
			}
		}

		// End of INIT
		
		// SETRHO - c  Converts position to computer normalization and accumulates charge density.
		std::vector<double> qdx(nsp);
		std::vector<double> rho(ng + 1, 0.0); // charge density at the spatial grid points
		std::vector<std::vector<double>> rhos(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> rho0(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> qjdata(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> qjp1data(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> drhodata(nsp, std::vector<double>(ng + 1, 0.0));

		double dx = gridx[1]; // spatial grid bin width
		double dtdx = dt / dx;

		for (int species = 0; species < nsp; species++) {
			qdx[species] = q[species] / dx;
			setRho(species, ng, dx, N, qdx, rho, x, rho0, rhos, iw);

			for (int K = 0; K < N[species]; ++K) {
				vx[species][K] *= dtdx;
			}
		}

		int t = 0;
		double ael = 1;

		// Acceleration
		std::vector<double> a(ng + 1);
		fields(rho, L, iw, dx, E, t, ng, a, ael);
		accel(nsp, dx, dt, t, q, m, ael, a, ng, N, x, vx);

		for (int i = 0; i < ng + 1; i++) {
			ESE[t] += std::pow(E[t][i],2) * 0.5 * dx;
		}

		Frame frame0;
		std::vector<Particle> particles;

		for (int species = 0; species < nsp; species++) {
			for (int i = 0; i < N[species]; i++) {
				int particleId = 0;

				Particle particle;

				particle.id = particleId;
				particle.position = x[species][i];
				particle.velocity = vx[species][i];
				particle.species = species;

				particles.push_back(particle);

				particleId++;
			}
		}

		frame0.particles = particles;
		frame0.electricField = E[t];
		frame0.frameNumber = t;
		mPicData.frames.push_back(frame0);

		//BEGIN TIME LOOP 

		nlohmann::json JSON;
		nlohmann::json frames;

		for (int t = 1; t <= nt; t++) {
			Frame frame;

			accel(nsp, dx, dt, t, q, m, ael, a, ng, N, x, vx);
			move(nsp, rho, rho0, qdx, N, x, vx, ng);
			fields(rho, L, iw, dx, E, t, ng, a, ael);

			for (int i = 0; i < ng + 1; i++) {
				ESE[t] += std::pow(E[t][i], 2) * 0.5 * dx;
			}

			frame.electricField = E[t];
			std::vector<Particle> particles;

			nlohmann::json JSONFrame;
			nlohmann::json JSONParticles;
			int particleId = 0;

			for (int species = 0; species < nsp; species++) {
				for (int i = 0; i < N[species]; i++) {

					Particle particle;


					particle.id = particleId;
					particle.position = x[species][i];
					particle.velocity = vx[species][i];
					ke[species][t] += 0.5*std::pow((particle.velocity), 2) * m[species];
					particle.species = species;

					particles.push_back(particle);


					nlohmann::json particleObject;

					particleObject["position"] = x[species][i];
					particleObject["velocity"] = vx[species][i];
					particleObject["species"] = species;
					particleObject["id"] = particleId;

					JSONParticles.push_back(particleObject);

					particleId++;
				}
			}

			JSONFrame["particles"] = JSONParticles;
			JSONFrame["frameNumber"] = t;
			frames.push_back(JSONFrame);

			frame.particles = particles;
			frame.frameNumber = t;
			mPicData.frames.push_back(frame);

			//// Create a scatter plot 
			//matplot::scatter(x[0], vx[0], 1);
			//matplot::hold(true);
			//matplot::scatter(x[1], vx[1], 1);

			//matplot::show();
			//matplot::hold(false);
			//std::cout << "t = " << t << std::endl;
		}

		JSON["ke"] = ke;
		JSON["ese"] = ESE;
		JSON["phaseFrames"] = frames;

		std::cout << JSON.dump() << std::endl;

		return true;
	}
};