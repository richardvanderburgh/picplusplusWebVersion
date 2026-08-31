#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include <DataStructs.h>
#include "Grid3D.hpp"

inline void depositParticle3d(
	const Grid3D& grid,
	std::vector<double>& rho,
	double charge,
	double x,
	double y,
	double z) {

	const double fx = x - std::floor(x);
	const double fy = y - std::floor(y);
	const double fz = z - std::floor(z);

	const int jx = static_cast<int>(std::floor(x));
	const int jy = static_cast<int>(std::floor(y));
	const int jz = static_cast<int>(std::floor(z));

	const int ix[2] = {grid.periodicIndex(jx, grid.nx), grid.periodicIndex(jx + 1, grid.nx)};
	const int iy[2] = {grid.periodicIndex(jy, grid.ny), grid.periodicIndex(jy + 1, grid.ny)};
	const int iz[2] = {grid.periodicIndex(jz, grid.nz), grid.periodicIndex(jz + 1, grid.nz)};

	const double wx[2] = {1.0 - fx, fx};
	const double wy[2] = {1.0 - fy, fy};
	const double wz[2] = {1.0 - fz, fz};

	for (int a = 0; a < 2; ++a) {
		for (int b = 0; b < 2; ++b) {
			for (int c = 0; c < 2; ++c) {
				const double weight = charge * wx[a] * wy[b] * wz[c];
				rho[grid.index(ix[a], iy[b], iz[c])] += weight;
			}
		}
	}
}

inline void setRho3d(
	int species,
	const DATA_STRUCTS::SimulationParams& simulationParams,
	const Grid3D& grid,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	std::vector<double>& rho) {

	if (species == 0) {
		std::fill(rho.begin(), rho.end(), 0.0);
	}

	const double invDx = 1.0 / grid.dx;
	const double invDy = 1.0 / grid.dy;
	const double invDz = 1.0 / grid.dz;
	// Match 1D qdx = q/dx: deposit charge density q / dV, not raw charge.
	const double chargeDensity = allSpeciesData[species].particleCharge / grid.cellVolume();

	auto& speciesData = allSpeciesData[species];
	for (int i = 0; i < speciesData.numParticles; ++i) {
		speciesData.particlePositions[i] *= invDx;
		speciesData.particlePositionsY[i] *= invDy;
		speciesData.particlePositionsZ[i] *= invDz;
		grid.wrapPosition(
			speciesData.particlePositions[i],
			speciesData.particlePositionsY[i],
			speciesData.particlePositionsZ[i]);
		depositParticle3d(
			grid,
			rho,
			chargeDensity,
			speciesData.particlePositions[i],
			speciesData.particlePositionsY[i],
			speciesData.particlePositionsZ[i]);
	}
}

inline void move3d(
	const DATA_STRUCTS::SimulationParams& simulationParams,
	const Grid3D& grid,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	std::vector<double>& rho) {

	std::fill(rho.begin(), rho.end(), 0.0);

	for (int species = 0; species < simulationParams.numSpecies; ++species) {
		std::vector<double>& px = allSpeciesData[species].particlePositions;
		std::vector<double>& py = allSpeciesData[species].particlePositionsY;
		std::vector<double>& pz = allSpeciesData[species].particlePositionsZ;
		const std::vector<double>& vx = allSpeciesData[species].particleXVelocities;
		const std::vector<double>& vy = allSpeciesData[species].particleYVelocities;
		const std::vector<double>& vz = allSpeciesData[species].particleZVelocities;
		const int numParticles = allSpeciesData[species].numParticles;
		const double chargeDensity = allSpeciesData[species].particleCharge / grid.cellVolume();

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			std::vector<double> localRho(rho.size(), 0.0);

#ifdef _OPENMP
#pragma omp for nowait schedule(static)
#endif
			for (int i = 0; i < numParticles; ++i) {
				double x = px[i] + vx[i];
				double y = py[i] + vy[i];
				double z = pz[i] + vz[i];
				grid.wrapPosition(x, y, z);
				px[i] = x;
				py[i] = y;
				pz[i] = z;

				const double fx = x - std::floor(x);
				const double fy = y - std::floor(y);
				const double fz = z - std::floor(z);
				const int jx = static_cast<int>(std::floor(x));
				const int jy = static_cast<int>(std::floor(y));
				const int jz = static_cast<int>(std::floor(z));
				const int ix[2] = {grid.periodicIndex(jx, grid.nx), grid.periodicIndex(jx + 1, grid.nx)};
				const int iy[2] = {grid.periodicIndex(jy, grid.ny), grid.periodicIndex(jy + 1, grid.ny)};
				const int iz[2] = {grid.periodicIndex(jz, grid.nz), grid.periodicIndex(jz + 1, grid.nz)};
				const double wx[2] = {1.0 - fx, fx};
				const double wy[2] = {1.0 - fy, fy};
				const double wz[2] = {1.0 - fz, fz};

				for (int a = 0; a < 2; ++a) {
					for (int b = 0; b < 2; ++b) {
						for (int c = 0; c < 2; ++c) {
							const double weight = chargeDensity * wx[a] * wy[b] * wz[c];
							localRho[grid.index(ix[a], iy[b], iz[c])] += weight;
						}
					}
				}
			}

#ifdef _OPENMP
#pragma omp critical
#endif
			{
				for (size_t k = 0; k < rho.size(); ++k) {
					rho[k] += localRho[k];
				}
			}
		}
	}
}
