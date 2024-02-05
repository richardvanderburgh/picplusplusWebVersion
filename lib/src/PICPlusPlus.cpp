
#include "PICPlusPlus.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <thread>
#include <vector>

#include "Accel.hpp"
#include "fft.hpp"
#include "Fields.hpp"
#include "SetRho.hpp"
#include "Utils.hpp"

namespace PIC_PLUS_PLUS {

	PICPlusPlus::PICPlusPlus(DATA_STRUCTS::InputVariables inputVariables)

		: m_simulationParams(inputVariables.simulationParams),
		m_allSpeciesData(inputVariables.allSpeciesData),

		m_qdx(inputVariables.simulationParams.numSpecies),

		m_rhos(inputVariables.simulationParams.numSpecies, std::vector<double>(m_simulationParams.numGrid + 1, 0.0)),
		m_rho0(m_simulationParams.numSpecies, std::vector<double>(m_simulationParams.numGrid + 1, 0.0)),

		m_chargeDensity(m_simulationParams.numGrid + 1, 0.0),
		m_electricField(m_simulationParams.numTimeSteps + 1, std::vector<double>(m_simulationParams.numGrid + 1, 0.0)),
		m_particleAcceleration(m_simulationParams.numGrid + 1),

		m_particleKineticEnergy(m_simulationParams.numSpecies, std::vector<double>(m_simulationParams.numTimeSteps + 1, 0.0)),
		m_electrostaticEnergy(m_simulationParams.numTimeSteps + 1, 0.0),
		m_totalEnergy(m_simulationParams.numTimeSteps + 1, 0.0)
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

			if (m_allSpeciesData[species].spatialPerturbationAmplitude != 0) {
				applySpatialPerturbation(m_allSpeciesData[species].particlePositions,
					m_allSpeciesData[species].numParticles,
					m_allSpeciesData[species].spatialPerturbationMode,
					m_allSpeciesData[species].spatialPerturbationAmplitude);
			}

			m_qdx[species] = m_allSpeciesData[species].particleCharge / m_simulationParams.gridStepSize;

			setRho(species,
				m_simulationParams,
				m_allSpeciesData,
				m_qdx,
				m_chargeDensity,
				m_rho0,
				m_rhos);
		}

		m_particleKineticEnergy.reserve(m_simulationParams.numSpecies);
		for (auto& speciesEnergy : m_particleKineticEnergy) {
			speciesEnergy.reserve(m_simulationParams.numTimeSteps + 1);
		}

		m_electrostaticEnergy.reserve(m_simulationParams.numTimeSteps + 1);
	};

	std::optional<nlohmann::json> PICPlusPlus::initialize() {

		fields(m_simulationParams,
			m_chargeDensity,
			m_electricField,
			m_timeStep,
			m_particleAcceleration,
			m_ael);

		accel(m_simulationParams,
			m_allSpeciesData,
			m_particleAcceleration,
			m_timeStep,
			m_ael);

		calculateEnergies();

		mPicData.frames.emplace_back() = updateFrame();

		auto start = std::chrono::high_resolution_clock::now();

		for (m_timeStep = 1; m_timeStep <= m_simulationParams.numTimeSteps; m_timeStep++) {
			runTimeLoop();
		}

		auto finish = std::chrono::high_resolution_clock::now();

		auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);

		std::cout << "Time loop took " << microseconds.count() << " micro secs\n";

		nlohmann::json JSON;
		JSON["ke"] = m_particleKineticEnergy;
		JSON["ese"] = m_electrostaticEnergy;
		JSON["phaseFrames"] = mPicData.frames;

		return JSON;
	}

	void PICPlusPlus::runTimeLoop() {

		accel(m_simulationParams,
			m_allSpeciesData,
			m_particleAcceleration,
			m_timeStep,
			m_ael);

		move(m_simulationParams,
			m_allSpeciesData,
			m_chargeDensity,
			m_rho0,
			m_qdx);

		fields(m_simulationParams,
			m_chargeDensity,
			m_electricField,
			m_timeStep,
			m_particleAcceleration,
			m_ael);

		calculateEnergies();

		mPicData.frames.emplace_back() = updateFrame();
	}

	void PICPlusPlus::initializePositions(std::vector<double>& inOutParticlePositions, const int numParticles, const double chargeCloudWidth) {
		for (int I = 1; I < numParticles + 1; I++) {
			inOutParticlePositions[I - 1] = (I - 0.5) * chargeCloudWidth;
		}
	}

	void PICPlusPlus::addDriftVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity) {
		for (int I = 0; I < numParticles; I++) {
			inOutParticleXVelocities[I] += driftVelocity;
		}
	}

	void PICPlusPlus::addThermalVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double thermalVelocity) {
		std::random_device rd;
		std::mt19937 gen(rd());
		std::normal_distribution<> dis(0, 1);

		for (int i = 0; i < numParticles; ++i) {
			double randomValue = dis(gen);
			inOutParticleXVelocities[i] += thermalVelocity * randomValue;
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

		double spatialPertConst = 2 * std::numbers::pi * spatialPerturbationMode / m_simulationParams.spatialLength;

		for (int a = 0; a < numParticles; ++a) {
			double theta = spatialPertConst * inOutParticlePositions[a];
			inOutParticlePositions[a] = inOutParticlePositions[a] + spatialPerturbationAmplitude * cos(theta);

			if (inOutParticlePositions[a] >= m_simulationParams.spatialLength) {
				inOutParticlePositions[a] = inOutParticlePositions[a] - m_simulationParams.spatialLength;
			}
			if (inOutParticlePositions[a] < 0) {
				inOutParticlePositions[a] = inOutParticlePositions[a] + m_simulationParams.spatialLength;
			}
		}
	}

	void PICPlusPlus::calculateEnergies() {

		for (int species = 0; species < m_simulationParams.numSpecies; species++) {
			for (int i = 0; i < m_allSpeciesData[species].numParticles; i++) {
				m_particleKineticEnergy[species][m_timeStep] += 0.5 * std::pow((m_allSpeciesData[species].particleXVelocities[i]), 2) * m_allSpeciesData[species].particleMass;
			}
		}

		double halfGridSize = 0.5 * m_simulationParams.gridStepSize;

		for (int i = 0; i < m_simulationParams.numGrid + 1; i++) {
			m_electrostaticEnergy[m_timeStep] += std::pow(m_electricField[m_timeStep][i], 2) * halfGridSize;
		}
	}

	DATA_STRUCTS::Frame PICPlusPlus::updateFrame() {

		DATA_STRUCTS::Frame frame;

		frame.electricField = m_electricField[m_timeStep];
		frame.particles = updateFrameParticles();
		frame.frameNumber = m_timeStep;

		return frame;
	}

	std::vector<DATA_STRUCTS::Particle> PICPlusPlus::updateFrameParticles() {

		int particleId = 0;

		std::vector<DATA_STRUCTS::Particle> particles;
		for (int species = 0; species < m_simulationParams.numSpecies; species++) {
			for (int i = 0; i < m_allSpeciesData[species].numParticles; i++) {

				DATA_STRUCTS::Particle& particle = particles.emplace_back();

				particle.id = particleId++;
				particle.position = m_allSpeciesData[species].particlePositions[i];
				particle.velocity = m_allSpeciesData[species].particleXVelocities[i];
				particle.species = species;
			}
		}
		return particles;
	}
};