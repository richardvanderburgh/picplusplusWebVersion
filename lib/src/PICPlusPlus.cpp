
#include "PICPlusPlus.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <thread>
#include <vector>

#include "Accel.hpp"
#include "PicPlusPlus.h"
#include "fft.hpp"
#include "Fields.hpp"
#include "SetRho.hpp"
#include "Utils.hpp"

namespace PIC_PLUS_PLUS {

	PICPlusPlus::PICPlusPlus(
		DATA_STRUCTS::SimulationParams simulationParams,
		std::vector<DATA_STRUCTS::SpeciesData> allSpeciesData)

		: m_simulationParams(simulationParams),
		m_allSpeciesData(allSpeciesData),

		m_qdx(simulationParams.numSpecies),

		m_rhos(simulationParams.numSpecies, std::vector<double>(m_simulationParams.numGrid + 1, 0.0)),
		m_rho0(m_simulationParams.numSpecies, std::vector<double>(m_simulationParams.numGrid + 1, 0.0)),

		m_particleAcceleration(m_simulationParams.numGrid + 1),

		m_electricField(m_simulationParams.numTimeSteps + 1, std::vector<double>(m_simulationParams.numGrid + 1, 0.0)),
		m_particleKineticEnergy(m_simulationParams.numSpecies, std::vector<double>(m_simulationParams.numTimeSteps + 1, 0.0)),
		m_electrostaticEnergy(m_simulationParams.numTimeSteps + 1, 0.0),
		m_totalEnergy(m_simulationParams.numTimeSteps + 1, 0.0),
		m_chargeDensity(m_simulationParams.numGrid + 1, 0.0)
	{
		m_simulationParams.gridStepSize = m_simulationParams.spatialLength / m_simulationParams.numGrid;
		m_dtdx = m_simulationParams.timeStepSize / m_simulationParams.gridStepSize;
		m_timeStep = 0;

		m_ael = 1;

		for (int species = 0; species < m_simulationParams.numSpecies; species++) {

			DATA_STRUCTS::SpeciesData& speciesData = m_allSpeciesData[species];

			speciesData.particleCharge = m_simulationParams.spatialLength * speciesData.plasmaFrequency * speciesData.plasmaFrequency / (speciesData.numParticles * speciesData.chargeMassRatio);
			speciesData.particleMass = speciesData.particleCharge / speciesData.chargeMassRatio;
			speciesData.chargeCloudWidth = m_simulationParams.spatialLength / speciesData.numParticles;

			initializePositions(speciesData.particlePositions, speciesData.numParticles, speciesData.chargeCloudWidth);
			initializeVelocities(speciesData.particleXVelocities, speciesData.numParticles, speciesData.driftVelocity, speciesData.thermalVelocity);

			for (int K = 0; K < speciesData.numParticles; ++K) {
				speciesData.particleXVelocities[K] *= m_dtdx;
			}

			m_particlePositions.push_back(speciesData.particlePositions);
			m_particleXVelocities.push_back(speciesData.particleXVelocities);
			m_particleCharge.push_back(speciesData.particleCharge);
			m_particleMass.push_back(speciesData.particleMass);
			m_speciesNumParticles.push_back(allSpeciesData[species].numParticles);

			if (m_allSpeciesData[species].spatialPerturbationAmplitude != 0) {
				applySpatialPerturbation(m_particlePositions[species],
					m_allSpeciesData[species].numParticles,
					m_allSpeciesData[species].spatialPerturbationMode,
					m_allSpeciesData[species].spatialPerturbationAmplitude);
			}

			m_qdx[species] = m_allSpeciesData[species].particleCharge / m_simulationParams.gridStepSize;

			setRho(species, m_simulationParams.numGrid, m_simulationParams.gridStepSize, m_allSpeciesData[species].numParticles,
				m_qdx, m_chargeDensity, m_particlePositions, m_rho0, m_rhos);
		}
	};

	std::optional<nlohmann::json> PICPlusPlus::initialize() {

		fields(m_simulationParams,
			m_chargeDensity,
			m_electricField,
			m_timeStep,
			m_particleAcceleration,
			m_ael);

		accel(
			m_simulationParams,
			m_particleAcceleration,
			m_particleCharge,
			m_timeStep,
			m_particleMass,
			m_ael,
			m_speciesNumParticles, 
			m_particlePositions, 
			m_particleXVelocities);

		calculateEnergies(
			m_particleXVelocities,
			m_particleMass);

		mPicData.frames.emplace_back() = updateFrame(
			m_simulationParams.numSpecies,
			m_speciesNumParticles,
			m_particlePositions,
			m_particleXVelocities,
			m_particleMass);

		for (m_timeStep = 1; m_timeStep <= m_simulationParams.numTimeSteps; m_timeStep++) {

			runTimeLoop();
		}

		nlohmann::json JSON;
		JSON["ke"] = m_particleKineticEnergy;
		JSON["ese"] = m_electrostaticEnergy;
		JSON["phaseFrames"] = mPicData.frames;

		return JSON;
	}

	void PICPlusPlus::runTimeLoop() {

		accel(m_simulationParams,
			m_particleAcceleration,
			m_particleCharge,
			m_timeStep,
			m_particleMass,
			m_ael,
			m_speciesNumParticles,
			m_particlePositions,
			m_particleXVelocities);

		move(m_simulationParams,
			m_chargeDensity,
			m_rho0,
			m_qdx, m_speciesNumParticles,
			m_particlePositions,
			m_particleXVelocities);

		fields(m_simulationParams,
			m_chargeDensity,
			m_electricField,
			m_timeStep,
			m_particleAcceleration,
			m_ael);

		calculateEnergies(
			m_particleXVelocities,
			m_particleMass);

		mPicData.frames.emplace_back() = updateFrame(
			m_simulationParams.numSpecies,
			m_speciesNumParticles,
			m_particlePositions,
			m_particleXVelocities,
			m_particleMass);
	}

	void PICPlusPlus::initializePositions(std::vector<double>& inOutParticlePositions, const int numParticles, const double chargeCloudWidth) {
		for (int I = 1; I < numParticles + 1; I++) {
			inOutParticlePositions[I - 1] = (I - 0.5) * chargeCloudWidth;
		}
	}

	void PICPlusPlus::addDriftVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity) {
		for (int I = 0; I < numParticles; I++) {
			inOutParticleXVelocities[I] = driftVelocity;
		}
	}

	void PICPlusPlus::addThermalVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double thermalVelocity) {

		for (int I = 0; I < numParticles; ++I) {
			double rm = 0;
			for (int ith = 0; ith < 12; ++ith) {
				rm += (double)rand() / RAND_MAX;
			}
			rm -= 6;

			inOutParticleXVelocities[I] += thermalVelocity * rm;
		}
	}

	void PICPlusPlus::initializeVelocities(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity, const double thermalVelocity) {
		if (driftVelocity != 0) {
			addDriftVelocity(inOutParticleXVelocities, numParticles, driftVelocity);
		}

		if (thermalVelocity != 0) {
			addThermalVelocity(inOutParticleXVelocities, numParticles, thermalVelocity);
		}
	}

	void PICPlusPlus::applySpatialPerturbation(std::vector<double>& inOutParticlePositions, const int numParticles, const int spatialPerturbationMode, const double spatialPerturbationAmplitude) {
		for (int a = 0; a < numParticles; ++a) {
			double theta = 2 * std::numbers::pi * spatialPerturbationMode * inOutParticlePositions[a] / m_simulationParams.spatialLength;
			inOutParticlePositions[a] = inOutParticlePositions[a] + spatialPerturbationAmplitude * cos(theta);

			if (inOutParticlePositions[a] >= m_simulationParams.spatialLength) {
				inOutParticlePositions[a] = inOutParticlePositions[a] - m_simulationParams.spatialLength;
			}
			if (inOutParticlePositions[a] < 0) {
				inOutParticlePositions[a] = inOutParticlePositions[a] + m_simulationParams.spatialLength;
			}
		}
	}

	void PICPlusPlus::calculateEnergies(
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass) {

		for (int species = 0; species < m_simulationParams.numSpecies; species++) {
			for (int i = 0; i < m_speciesNumParticles[species]; i++) {
				m_particleKineticEnergy[species][m_timeStep] += 0.5 * std::pow((particleXVelocities[species][i]), 2) * particleMass[species];
			}
		}

		for (int i = 0; i < m_simulationParams.numGrid + 1; i++) {
			m_electrostaticEnergy[m_timeStep] += std::pow(m_electricField[m_timeStep][i], 2) * 0.5 * m_simulationParams.gridStepSize;
		}
	}

	std::vector<DATA_STRUCTS::Particle> PICPlusPlus::updateFrameParticles(
		int numSpecies,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particlePositions,
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass) {

		int particleId = 0;

		std::vector<DATA_STRUCTS::Particle> particles;
		for (int species = 0; species < numSpecies; species++) {
			for (int i = 0; i < speciesNumParticles[species]; i++) {

				DATA_STRUCTS::Particle& particle = particles.emplace_back();

				particle.id = particleId++;
				particle.position = particlePositions[species][i];
				particle.velocity = particleXVelocities[species][i];
				particle.species = species;
			}
		}
		return particles;
	}

	DATA_STRUCTS::Frame PICPlusPlus::updateFrame(
		const int numSpecies,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particlePositions,
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass) {

		DATA_STRUCTS::Frame frame;

		frame.electricField = m_electricField[m_timeStep];

		frame.particles = updateFrameParticles(
			numSpecies,
			speciesNumParticles,
			particlePositions,
			particleXVelocities,
			particleMass
		);

		frame.frameNumber = m_timeStep;

		return frame;
	}
};