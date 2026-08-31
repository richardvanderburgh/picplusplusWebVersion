#pragma once

#include <vector>

#include "complex.hpp"
#include "fft.hpp"
#include "Grid3D.hpp"

// Separable 3D complex FFT using the bundled 1D radix-2 FFT.
// Data layout: index(i, j, k) = i + nx * (j + ny * k), i in [0, nx), etc.
inline bool fft3dForward(std::vector<complex>& data, const Grid3D& grid) {
	const int nx = grid.nx;
	const int ny = grid.ny;
	const int nz = grid.nz;
	if (static_cast<int>(data.size()) != nx * ny * nz) {
		return false;
	}

	std::vector<complex> line(std::max({nx, ny, nz}));

	for (int j = 0; j < ny; ++j) {
		for (int k = 0; k < nz; ++k) {
			for (int i = 0; i < nx; ++i) {
				line[static_cast<size_t>(i)] = data[grid.index(i, j, k)];
			}
			if (!CFFT::Forward(line.data(), static_cast<unsigned int>(nx))) {
				return false;
			}
			for (int i = 0; i < nx; ++i) {
				data[grid.index(i, j, k)] = line[static_cast<size_t>(i)];
			}
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 0; k < nz; ++k) {
			for (int j = 0; j < ny; ++j) {
				line[static_cast<size_t>(j)] = data[grid.index(i, j, k)];
			}
			if (!CFFT::Forward(line.data(), static_cast<unsigned int>(ny))) {
				return false;
			}
			for (int j = 0; j < ny; ++j) {
				data[grid.index(i, j, k)] = line[static_cast<size_t>(j)];
			}
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (int k = 0; k < nz; ++k) {
				line[static_cast<size_t>(k)] = data[grid.index(i, j, k)];
			}
			if (!CFFT::Forward(line.data(), static_cast<unsigned int>(nz))) {
				return false;
			}
			for (int k = 0; k < nz; ++k) {
				data[grid.index(i, j, k)] = line[static_cast<size_t>(k)];
			}
		}
	}

	return true;
}

inline bool fft3dInverse(std::vector<complex>& data, const Grid3D& grid) {
	const int nx = grid.nx;
	const int ny = grid.ny;
	const int nz = grid.nz;
	if (static_cast<int>(data.size()) != nx * ny * nz) {
		return false;
	}

	std::vector<complex> line(std::max({nx, ny, nz}));

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (int k = 0; k < nz; ++k) {
				line[static_cast<size_t>(k)] = data[grid.index(i, j, k)];
			}
			if (!CFFT::Inverse(line.data(), static_cast<unsigned int>(nz))) {
				return false;
			}
			for (int k = 0; k < nz; ++k) {
				data[grid.index(i, j, k)] = line[static_cast<size_t>(k)];
			}
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 0; k < nz; ++k) {
			for (int j = 0; j < ny; ++j) {
				line[static_cast<size_t>(j)] = data[grid.index(i, j, k)];
			}
			if (!CFFT::Inverse(line.data(), static_cast<unsigned int>(ny))) {
				return false;
			}
			for (int j = 0; j < ny; ++j) {
				data[grid.index(i, j, k)] = line[static_cast<size_t>(j)];
			}
		}
	}

	for (int j = 0; j < ny; ++j) {
		for (int k = 0; k < nz; ++k) {
			for (int i = 0; i < nx; ++i) {
				line[static_cast<size_t>(i)] = data[grid.index(i, j, k)];
			}
			if (!CFFT::Inverse(line.data(), static_cast<unsigned int>(nx))) {
				return false;
			}
			for (int i = 0; i < nx; ++i) {
				data[grid.index(i, j, k)] = line[static_cast<size_t>(i)];
			}
		}
	}

	return true;
}

inline int fftMode(int index, int half) {
	if (index == 0) {
		return 0;
	}
	if (index <= half) {
		return index;
	}
	return index - 2 * half;
}
