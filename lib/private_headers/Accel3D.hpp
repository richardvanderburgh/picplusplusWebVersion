#pragma once

#include <cmath>
#include <vector>

#include <DataStructs.h>
#include "Boris.hpp"
#include "Grid3D.hpp"

inline double trilinearSample(
	const std::vector<double>& field,
	const Grid3D& grid,
	double x,
	double y,
	double z) {

	const double fx = x - std::floor(x);
	const double fy = y - std::floor(y);
	const double fz = z - std::floor(z);

	const int jx = static_cast<int>(std::floor(x));
	const int jy = static_cast<int>(std::floor(y));
	const int jz = static_cast<int>(std::floor(z));

	const int jxp = grid.periodicIndex(jx + 1, grid.nx);
	const int jyp = grid.periodicIndex(jy + 1, grid.ny);
	const int jzp = grid.periodicIndex(jz + 1, grid.nz);
	const int jxm = grid.periodicIndex(jx, grid.nx);
	const int jym = grid.periodicIndex(jy, grid.ny);
	const int jzm = grid.periodicIndex(jz, grid.nz);

	const double c000 = field[grid.index(jxm, jym, jzm)];
	const double c100 = field[grid.index(jxp, jym, jzm)];
	const double c010 = field[grid.index(jxm, jyp, jzm)];
	const double c110 = field[grid.index(jxp, jyp, jzm)];
	const double c001 = field[grid.index(jxm, jym, jzp)];
	const double c101 = field[grid.index(jxp, jym, jzp)];
	const double c011 = field[grid.index(jxm, jyp, jzp)];
	const double c111 = field[grid.index(jxp, jyp, jzp)];

	const double c00 = c000 * (1.0 - fx) + c100 * fx;
	const double c10 = c010 * (1.0 - fx) + c110 * fx;
	const double c01 = c001 * (1.0 - fx) + c101 * fx;
	const double c11 = c011 * (1.0 - fx) + c111 * fx;
	const double c0 = c00 * (1.0 - fy) + c10 * fy;
	const double c1 = c01 * (1.0 - fy) + c11 * fy;
	return c0 * (1.0 - fz) + c1 * fz;
}

inline void accel3d(
	const DATA_STRUCTS::SimulationParams& simulationParams,
	const Grid3D& grid,
	std::vector<DATA_STRUCTS::SpeciesData>& allSpeciesData,
	const std::vector<double>& ex,
	const std::vector<double>& ey,
	const std::vector<double>& ez,
	int timeStep,
	double& ael) {

	const double dxdt = grid.dx / simulationParams.timeStepSize;
	const double dydt = grid.dy / simulationParams.timeStepSize;
	const double dzdt = grid.dz / simulationParams.timeStepSize;
	const bool magnetized = simulationParams.hasMagneticField();
	const double bx = simulationParams.magneticFieldX;
	const double by = simulationParams.magneticFieldY;
	const double bz = simulationParams.magneticFieldZ;
	const double dt = simulationParams.timeStepSize;
	// Leapfrog half-step back at t = 0 (electric + magnetic).
	const double pushDt = (timeStep == 0) ? -0.5 * dt : dt;

	for (int species = 0; species < simulationParams.numSpecies; ++species) {
		const double qm = allSpeciesData[species].particleCharge / allSpeciesData[species].particleMass;
		// ae converts E → full (or half at t=0) code-unit Δv for this push.
		const double aeX = qm * (pushDt / dxdt);
		const double aeY = qm * (pushDt / dydt);
		const double aeZ = qm * (pushDt / dzdt);
		ael = aeX;

		std::vector<double>& vx = allSpeciesData[species].particleXVelocities;
		std::vector<double>& vy = allSpeciesData[species].particleYVelocities;
		std::vector<double>& vz = allSpeciesData[species].particleZVelocities;
		const std::vector<double>& px = allSpeciesData[species].particlePositions;
		const std::vector<double>& py = allSpeciesData[species].particlePositionsY;
		const std::vector<double>& pz = allSpeciesData[species].particlePositionsZ;
		const int numParticles = allSpeciesData[species].numParticles;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
		for (int i = 0; i < numParticles; ++i) {
			const double ax = trilinearSample(ex, grid, px[i], py[i], pz[i]);
			const double ay = trilinearSample(ey, grid, px[i], py[i], pz[i]);
			const double az = trilinearSample(ez, grid, px[i], py[i], pz[i]);

			if (!magnetized) {
				vx[i] += aeX * ax;
				vy[i] += aeY * ay;
				vz[i] += aeZ * az;
				continue;
			}

			borisPushCodeUnits(
				vx[i], vy[i], vz[i],
				ax, ay, az,
				aeX, aeY, aeZ,
				qm, pushDt,
				dxdt, dydt, dzdt,
				bx, by, bz);
		}
	}
}
