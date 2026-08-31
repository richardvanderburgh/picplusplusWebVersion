#pragma once

#include <cmath>
#include <cstdint>
#include <numbers>
#include <stdexcept>
#include <vector>

#include <DataStructs.h>
#include "FFT3D.hpp"
#include "Grid3D.hpp"

inline void fields3d(
	const Grid3D& grid,
	std::vector<double>& rho,
	std::vector<double>& ex,
	std::vector<double>& ey,
	std::vector<double>& ez) {

	const int nx = grid.nx;
	const int ny = grid.ny;
	const int nz = grid.nz;
	const size_t nCells = grid.numCells();

	if (rho.size() != nCells || ex.size() != nCells || ey.size() != nCells || ez.size() != nCells) {
		throw std::runtime_error("fields3d: grid buffer size mismatch");
	}

	std::vector<complex> rhoK(nCells);
	std::vector<complex> phiK(nCells);
	std::vector<double> electricPotential(nCells);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nCells >= 4096)
#endif
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				rhoK[grid.index(i, j, k)] = rho[grid.index(i, j, k)];
			}
		}
	}

	if (!fft3dForward(rhoK, grid)) {
		throw std::runtime_error("fields3d: forward FFT failed");
	}

	const double twoPi = 2.0 * std::numbers::pi;
	const int halfX = nx / 2;
	const int halfY = ny / 2;
	const int halfZ = nz / 2;

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nCells >= 4096)
#endif
	for (int k = 0; k < nz; ++k) {
		const int mz = fftMode(k, halfZ);
		const double kz = twoPi * static_cast<double>(mz) / grid.lz;
		for (int j = 0; j < ny; ++j) {
			const int my = fftMode(j, halfY);
			const double ky = twoPi * static_cast<double>(my) / grid.ly;
			for (int i = 0; i < nx; ++i) {
				const int mx = fftMode(i, halfX);
				const double kx = twoPi * static_cast<double>(mx) / grid.lx;
				const double kSquared = kx * kx + ky * ky + kz * kz;
				const size_t idx = grid.index(i, j, k);
				if (kSquared <= 0.0) {
					phiK[idx] = complex(0.0, 0.0);
				} else {
					phiK[idx] = rhoK[idx] / (-kSquared);
				}
			}
		}
	}

	if (!fft3dInverse(phiK, grid)) {
		throw std::runtime_error("fields3d: inverse FFT failed");
	}

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nCells >= 4096)
#endif
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				electricPotential[grid.index(i, j, k)] = -phiK[grid.index(i, j, k)].re();
			}
		}
	}

	const double invTwoDx = 1.0 / (2.0 * grid.dx);
	const double invTwoDy = 1.0 / (2.0 * grid.dy);
	const double invTwoDz = 1.0 / (2.0 * grid.dz);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) if(nCells >= 4096)
#endif
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				const int ip = grid.periodicIndex(i + 1, nx);
				const int im = grid.periodicIndex(i - 1, nx);
				const int jp = grid.periodicIndex(j + 1, ny);
				const int jm = grid.periodicIndex(j - 1, ny);
				const int kp = grid.periodicIndex(k + 1, nz);
				const int km = grid.periodicIndex(k - 1, nz);
				const size_t idx = grid.index(i, j, k);

				ex[idx] = (electricPotential[grid.index(im, j, k)] - electricPotential[grid.index(ip, j, k)]) * invTwoDx;
				ey[idx] = (electricPotential[grid.index(i, jm, k)] - electricPotential[grid.index(i, jp, k)]) * invTwoDy;
				ez[idx] = (electricPotential[grid.index(i, j, km)] - electricPotential[grid.index(i, j, kp)]) * invTwoDz;
			}
		}
	}
}
