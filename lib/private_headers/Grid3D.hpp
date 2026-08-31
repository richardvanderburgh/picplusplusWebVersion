#pragma once

#include <cmath>
#include <cstddef>

#include <DataStructs.h>

struct Grid3D {
	int nx = 0;
	int ny = 0;
	int nz = 0;
	double lx = 0.0;
	double ly = 0.0;
	double lz = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;

	static Grid3D fromParams(const DATA_STRUCTS::SimulationParams& params) {
		Grid3D grid;
		grid.nx = params.numGrid;
		grid.ny = params.numGridY > 0 ? params.numGridY : params.numGrid;
		grid.nz = params.numGridZ > 0 ? params.numGridZ : params.numGrid;
		grid.lx = params.spatialLength;
		grid.ly = params.spatialLengthY > 0.0 ? params.spatialLengthY : params.spatialLength;
		grid.lz = params.spatialLengthZ > 0.0 ? params.spatialLengthZ : params.spatialLength;
		grid.dx = grid.lx / grid.nx;
		grid.dy = grid.ly / grid.ny;
		grid.dz = grid.lz / grid.nz;
		return grid;
	}

	[[nodiscard]] size_t numCells() const {
		return static_cast<size_t>(nx) * static_cast<size_t>(ny) * static_cast<size_t>(nz);
	}

	[[nodiscard]] size_t index(int i, int j, int k) const {
		return static_cast<size_t>(i)
			+ static_cast<size_t>(nx) * (static_cast<size_t>(j) + static_cast<size_t>(ny) * static_cast<size_t>(k));
	}

	[[nodiscard]] double cellVolume() const {
		return dx * dy * dz;
	}

	static void wrapGridPosition(double& position, int numCells) {
		const double length = static_cast<double>(numCells);
		while (position < 0.0) {
			position += length;
		}
		while (position >= length) {
			position -= length;
		}
	}

	void wrapPosition(double& x, double& y, double& z) const {
		wrapGridPosition(x, nx);
		wrapGridPosition(y, ny);
		wrapGridPosition(z, nz);
	}

	[[nodiscard]] int periodicIndex(int i, int axisSize) const {
		int wrapped = i % axisSize;
		if (wrapped < 0) {
			wrapped += axisSize;
		}
		return wrapped;
	}
};
