#include "PICPlusPlus.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <random>

#include "Accel3D.hpp"
#include "Fields3D.hpp"
#include "SetRho3D.hpp"
#include "Utils.hpp"

namespace PIC_PLUS_PLUS {

namespace {

void latticeCounts(int numParticles, int& nx, int& ny, int& nz) {
	nx = std::max(1, static_cast<int>(std::cbrt(static_cast<double>(numParticles))));
	ny = std::max(1, static_cast<int>(std::sqrt(static_cast<double>(numParticles) / nx)));
	nz = std::max(1, numParticles / (nx * ny));
	while (nx * ny * nz < numParticles) {
		++nz;
	}
}

} // namespace

void PICPlusPlus::initializePositions3D(DATA_STRUCTS::SpeciesData& speciesData) {
	int nx = 0;
	int ny = 0;
	int nz = 0;
	latticeCounts(speciesData.numParticles, nx, ny, nz);

	const double widthX = m_grid3D.lx / nx;
	const double widthY = m_grid3D.ly / ny;
	const double widthZ = m_grid3D.lz / nz;

	for (int p = 0; p < speciesData.numParticles; ++p) {
		const int iz = p / (nx * ny);
		const int iy = (p % (nx * ny)) / nx;
		const int ix = p % nx;
		speciesData.particlePositions[static_cast<size_t>(p)] = (ix + 0.5) * widthX;
		speciesData.particlePositionsY[static_cast<size_t>(p)] = (iy + 0.5) * widthY;
		speciesData.particlePositionsZ[static_cast<size_t>(p)] = (iz + 0.5) * widthZ;
	}
}

void PICPlusPlus::addThermalVelocity3D(DATA_STRUCTS::SpeciesData& speciesData) {
	std::mt19937 gen(42);
	std::normal_distribution<> dis(0, 1);
	const double vth = speciesData.thermalVelocity;

	for (int i = 0; i < speciesData.numParticles; ++i) {
		speciesData.particleXVelocities[static_cast<size_t>(i)] += vth * dis(gen);
		speciesData.particleYVelocities[static_cast<size_t>(i)] += vth * dis(gen);
		speciesData.particleZVelocities[static_cast<size_t>(i)] += vth * dis(gen);
	}
}

void PICPlusPlus::initializeVelocities3D(DATA_STRUCTS::SpeciesData& speciesData) {
	if (speciesData.driftVelocity != 0.0) {
		for (int i = 0; i < speciesData.numParticles; ++i) {
			speciesData.particleXVelocities[static_cast<size_t>(i)] += speciesData.driftVelocity;
		}
	}
	if (speciesData.thermalVelocity != 0.0) {
		addThermalVelocity3D(speciesData);
	}

	const double dtdx = m_simulationParams.timeStepSize / m_grid3D.dx;
	const double dtdy = m_simulationParams.timeStepSize / m_grid3D.dy;
	const double dtdz = m_simulationParams.timeStepSize / m_grid3D.dz;
	for (int i = 0; i < speciesData.numParticles; ++i) {
		speciesData.particleXVelocities[static_cast<size_t>(i)] *= dtdx;
		speciesData.particleYVelocities[static_cast<size_t>(i)] *= dtdy;
		speciesData.particleZVelocities[static_cast<size_t>(i)] *= dtdz;
	}
}

void PICPlusPlus::applySpatialPerturbation3D(DATA_STRUCTS::SpeciesData& speciesData) {
	const double spatialPertConst = 2.0 * std::numbers::pi * speciesData.spatialPerturbationMode / m_grid3D.lx;
	const bool useSin = speciesData.spatialPerturbationWaveform == "sin";

	for (int a = 0; a < speciesData.numParticles; ++a) {
		const double theta = spatialPertConst * speciesData.particlePositions[static_cast<size_t>(a)];
		const double shift = speciesData.spatialPerturbationAmplitude
			* (useSin ? std::sin(theta) : std::cos(theta));
		speciesData.particlePositions[static_cast<size_t>(a)] += shift;

		if (speciesData.particlePositions[static_cast<size_t>(a)] >= m_grid3D.lx) {
			speciesData.particlePositions[static_cast<size_t>(a)] -= m_grid3D.lx;
		}
		if (speciesData.particlePositions[static_cast<size_t>(a)] < 0.0) {
			speciesData.particlePositions[static_cast<size_t>(a)] += m_grid3D.lx;
		}
	}
}

void PICPlusPlus::calculateEnergies3D() {
	const double cellVolume = m_grid3D.cellVolume();
	// Velocities are stored in cells/timestep; convert to physical units so
	// KE is commensurate with PE = (1/2) ∫ E² dV.
	const double dxdt = m_grid3D.dx / m_simulationParams.timeStepSize;
	const double dydt = m_grid3D.dy / m_simulationParams.timeStepSize;
	const double dzdt = m_grid3D.dz / m_simulationParams.timeStepSize;

	for (int species = 0; species < m_simulationParams.numSpecies; ++species) {
		const double particleMass = m_allSpeciesData[species].particleMass;
		const int numParticles = m_allSpeciesData[species].numParticles;
		double kineticEnergy = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : kineticEnergy) schedule(static)
#endif
		for (int i = 0; i < numParticles; ++i) {
			const double vx = m_allSpeciesData[species].particleXVelocities[static_cast<size_t>(i)] * dxdt;
			const double vy = m_allSpeciesData[species].particleYVelocities[static_cast<size_t>(i)] * dydt;
			const double vz = m_allSpeciesData[species].particleZVelocities[static_cast<size_t>(i)] * dzdt;
			kineticEnergy += 0.5 * (vx * vx + vy * vy + vz * vz) * particleMass;
		}
		m_particleKineticEnergy[species][m_timeStep] += kineticEnergy;
	}

	const size_t nCells = m_grid3D.numCells();
	double instantaneousFieldEnergy = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : instantaneousFieldEnergy) schedule(static) if(nCells >= 4096)
#endif
	for (size_t idx = 0; idx < nCells; ++idx) {
		const double ex = m_electricFieldX3D[idx];
		const double ey = m_electricFieldY3D[idx];
		const double ez = m_electricFieldZ3D[idx];
		instantaneousFieldEnergy += (ex * ex + ey * ey + ez * ez) * 0.5 * cellVolume;
	}
	// Pair KE(v^{n+1/2}) with ½(W_E^n + W_E^{n+1}) to remove the artificial
	// 2ω_p wobble from mixing half-step KE with a single integer-step PE.
	const double reportedFieldEnergy = (m_timeStep == 0)
		? instantaneousFieldEnergy
		: 0.5 * (m_previousFieldEnergy + instantaneousFieldEnergy);
	m_electrostaticEnergy[m_timeStep] += reportedFieldEnergy;
	m_previousFieldEnergy = instantaneousFieldEnergy;
}

std::vector<DATA_STRUCTS::Particle> PICPlusPlus::updateFrameParticles3D() {
	std::vector<DATA_STRUCTS::Particle> particles;
	int particleId = 0;

	const double dtdx = m_simulationParams.timeStepSize / m_grid3D.dx;
	const double dtdy = m_simulationParams.timeStepSize / m_grid3D.dy;
	const double dtdz = m_simulationParams.timeStepSize / m_grid3D.dz;

	for (int species = 0; species < m_simulationParams.numSpecies; ++species) {
		for (int i = 0; i < m_allSpeciesData[species].numParticles; ++i) {
			DATA_STRUCTS::Particle& particle = particles.emplace_back();
			particle.id = particleId++;
			particle.species = species;
			particle.position = m_allSpeciesData[species].particlePositions[static_cast<size_t>(i)] * m_grid3D.dx;
			particle.positionY = m_allSpeciesData[species].particlePositionsY[static_cast<size_t>(i)] * m_grid3D.dy;
			particle.positionZ = m_allSpeciesData[species].particlePositionsZ[static_cast<size_t>(i)] * m_grid3D.dz;
			particle.velocity = m_allSpeciesData[species].particleXVelocities[static_cast<size_t>(i)] / dtdx;
			particle.velocityY = m_allSpeciesData[species].particleYVelocities[static_cast<size_t>(i)] / dtdy;
			particle.velocityZ = m_allSpeciesData[species].particleZVelocities[static_cast<size_t>(i)] / dtdz;
		}
	}
	return particles;
}

DATA_STRUCTS::Frame PICPlusPlus::updateFrame3D() {
	DATA_STRUCTS::Frame frame;
	frame.dimension = 3;
	frame.particles = updateFrameParticles3D();
	frame.frameNumber = m_timeStep;

	const int midY = m_grid3D.ny / 2;
	const int midZ = m_grid3D.nz / 2;
	frame.electricField.resize(static_cast<size_t>(m_grid3D.nx));
	for (int i = 0; i < m_grid3D.nx; ++i) {
		frame.electricField[static_cast<size_t>(i)] =
			m_electricFieldX3D[m_grid3D.index(i, midY, midZ)];
	}
	return frame;
}

void PICPlusPlus::runTimeLoop3D() {
	accel3d(
		m_simulationParams,
		m_grid3D,
		m_allSpeciesData,
		m_electricFieldX3D,
		m_electricFieldY3D,
		m_electricFieldZ3D,
		m_timeStep,
		m_ael);

	move3d(m_simulationParams, m_grid3D, m_allSpeciesData, m_chargeDensity3D);

	fields3d(
		m_grid3D,
		m_chargeDensity3D,
		m_electricFieldX3D,
		m_electricFieldY3D,
		m_electricFieldZ3D);

	calculateEnergies3D();

	if (m_simulationParams.framePeriod > 0 && m_timeStep % m_simulationParams.framePeriod == 0) {
		mPicData.frames.emplace_back() = updateFrame3D();
	}
}

std::optional<nlohmann::json> PICPlusPlus::initialize3D() {
	if (const auto validationError = validateSimulationParams(m_simulationParams)) {
		std::cerr << "Invalid simulation parameters: " << *validationError << "\n";
		return std::nullopt;
	}

	fields3d(
		m_grid3D,
		m_chargeDensity3D,
		m_electricFieldX3D,
		m_electricFieldY3D,
		m_electricFieldZ3D);

	accel3d(
		m_simulationParams,
		m_grid3D,
		m_allSpeciesData,
		m_electricFieldX3D,
		m_electricFieldY3D,
		m_electricFieldZ3D,
		m_timeStep,
		m_ael);

	calculateEnergies3D();
	mPicData.frames.emplace_back() = updateFrame3D();

	const int totalSteps = m_simulationParams.numTimeSteps;
	const bool reportProgress = std::getenv("PICPP_PROGRESS") != nullptr;
	int lastReportedPercent = -1;

	auto start = std::chrono::high_resolution_clock::now();

	for (m_timeStep = 1; m_timeStep <= totalSteps; m_timeStep++) {
		runTimeLoop3D();

		if (reportProgress && totalSteps > 0) {
			const int percent = (100 * m_timeStep) / totalSteps;
			if (percent != lastReportedPercent) {
				std::cerr << "PICPP_PROGRESS " << percent << '\n' << std::flush;
				lastReportedPercent = percent;
			}
		}
	}

	if (reportProgress) {
		std::cerr << "PICPP_PROGRESS 100\n" << std::flush;
	}

	auto finish = std::chrono::high_resolution_clock::now();
	auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
	std::cout << "Time loop took " << microseconds.count() << " micro secs\n";

	nlohmann::json jsonResult;
	jsonResult["dimension"] = 3;
	jsonResult["ke"] = m_particleKineticEnergy;
	jsonResult["ese"] = m_electrostaticEnergy;
	jsonResult["phaseFrames"] = mPicData.frames;
	return jsonResult;
}

} // namespace PIC_PLUS_PLUS
