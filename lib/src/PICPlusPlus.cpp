
#include "PICPlusPlus.h"

//#include <fftw3.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <thread>
#include <vector>

#include "Accel.hpp"
#include "DataStructs.h"
#include "fft.hpp"
#include "Fields.hpp"
#include "SetRho.hpp"
#include "Utils.hpp"


namespace PIC_PLUS_PLUS {

	PICPlusPlus::PICPlusPlus(
		const double spatialLength,
		const int numParticles,
		const int numTimeSteps,
		const double timeStepSize,
		const int numGrid,
		const int spatialPerturbationMode,
		const double driftVelocity,
		const int numSpecies,
		const double spatialPerturbationAmplitude,
		const double thermalVelocity,
		const int plasmaFrequency,
		const int chargeMassRatio)

		: m_spatialLength(spatialLength),
		m_numParticles(numParticles),
		m_numTimeSteps(numTimeSteps),
		m_timeStepSize(timeStepSize),
		m_numGrid(numGrid),

		m_spatialPerturbationMode(spatialPerturbationMode),
		m_driftVelocity(driftVelocity),
		m_numSpecies(numSpecies),
		m_spatialPerturbationAmplitude(spatialPerturbationAmplitude),
		m_thermalVelocity(thermalVelocity),
		m_plasmaFrequency(plasmaFrequency),
		m_chargeMassRatio(chargeMassRatio),

	    m_chargeCloudWidth(m_numSpecies),
	    m_particleCharge(m_numSpecies, 0.0),
	    m_particleMass(m_numSpecies),
	    m_qdx(m_numSpecies),

		m_speciesNumParticles(m_numSpecies, m_numParticles),

		m_rhos(m_numSpecies, std::vector<double>(m_numGrid + 1, 0.0)),
		m_rho0(m_numSpecies, std::vector<double>(m_numGrid + 1, 0.0)),


		m_particleAcceleration(m_numGrid + 1),

		m_electricField(m_numTimeSteps + 1, std::vector<double>(m_numGrid + 1, 0.0)),
		m_particleKineticEnergy(m_numSpecies, std::vector<double>(m_numTimeSteps + 1, 0.0)),
		m_electrostaticEnergy(m_numTimeSteps + 1, 0.0),
		m_totalEnergy(m_numTimeSteps + 1, 0.0),
		m_chargeDensity(m_numGrid + 1, 0.0),
	
		m_particlePositions(m_numSpecies, std::vector<double>(m_numSpecies * m_numParticles, 0.0)),
		m_particleXVelocities(m_numSpecies, std::vector<double>(m_numSpecies * m_numParticles, 0.0))
	{

		m_gridStepSize = m_spatialLength / m_numGrid;
		m_dtdx = m_timeStepSize / m_gridStepSize;
		m_timeStep = 0;

		m_ael = 1;

	};

	std::optional<nlohmann::json> PICPlusPlus::initialize() {
 
		std::vector<int>		  speciesSpatialPerturbationMode(m_numSpecies, m_spatialPerturbationMode);
		std::vector<double>		  speciesSpatialPerturbationAmplitude(m_numSpecies, m_spatialPerturbationAmplitude);

		std::vector<double> speciesPlasmaFrequency(m_numSpecies, m_plasmaFrequency);
		std::vector<double> speciesChargeMassRatio(m_numSpecies, m_chargeMassRatio);

		const std::vector<double> speciesDriftVelocity = { m_driftVelocity, - m_driftVelocity, 0, 2 *  m_driftVelocity, -2 *  m_driftVelocity };

		std::vector<double> speciesThermalVelocity(m_numSpecies, m_thermalVelocity);

		for (int species = 0; species < m_numSpecies; species++) {

			initializeSpecies(m_particleCharge[species], 
				m_particleMass[species], 
				m_chargeCloudWidth[species], 
				m_particlePositions[species], 
				m_particleXVelocities[species],
				speciesPlasmaFrequency[species],
				m_speciesNumParticles[species],
				speciesChargeMassRatio[species],
				speciesPlasmaFrequency[species],
				speciesDriftVelocity[species],
				speciesThermalVelocity[species]);

			if (speciesSpatialPerturbationAmplitude[species] != 0) {
				applySpatialPerturbation(m_particlePositions[species], 
					m_speciesNumParticles[species], 
					speciesSpatialPerturbationMode[species],
					speciesSpatialPerturbationAmplitude[species]);
			}

			m_qdx[species] = m_particleCharge[species] / m_gridStepSize;
			setRho(species, m_numGrid, m_gridStepSize, m_speciesNumParticles, m_qdx, m_chargeDensity, m_particlePositions, m_rho0, m_rhos);
		}

		fields(m_chargeDensity, m_spatialLength, m_gridStepSize, m_electricField, m_timeStep, m_numGrid, m_particleAcceleration, m_ael);

		accel(m_particleAcceleration,
			m_particleCharge,
			m_numSpecies,
			m_gridStepSize,
			m_timeStepSize,
			m_timeStep,
			m_particleMass,
			m_ael, m_numGrid, m_speciesNumParticles, m_particlePositions, m_particleXVelocities);

		calculateEnergies(
			m_particleXVelocities,
			m_particleMass);

		mPicData.frames.emplace_back() = updateFrame(
			m_numSpecies,
			m_speciesNumParticles,
			m_particlePositions,
			m_particleXVelocities,
			m_particleMass);

		//end of initializing
		
		for (m_timeStep = 1; m_timeStep <= m_numTimeSteps; m_timeStep++) {

			runTimeLoop();
		}

		nlohmann::json JSON;
		JSON["ke"] = m_particleKineticEnergy;
		JSON["ese"] = m_electrostaticEnergy;
		JSON["phaseFrames"] = mPicData.frames;

		return JSON;
	}

	void PICPlusPlus::runTimeLoop() {
		
			accel(m_particleAcceleration,
				m_particleCharge,
				m_numSpecies,
				m_gridStepSize,
				m_timeStepSize,
				m_timeStep,
				m_particleMass,
				m_ael, m_numGrid, m_speciesNumParticles, m_particlePositions, m_particleXVelocities);

			move(m_numSpecies, m_chargeDensity, m_rho0, m_qdx, m_speciesNumParticles, m_particlePositions, m_particleXVelocities, m_numGrid);
			fields(m_chargeDensity, m_spatialLength, m_gridStepSize, m_electricField, m_timeStep, m_numGrid, m_particleAcceleration, m_ael);

			calculateEnergies(
				m_particleXVelocities,
				m_particleMass);

			mPicData.frames.emplace_back() = updateFrame(
				m_numSpecies,
				m_speciesNumParticles,
				m_particlePositions,
				m_particleXVelocities,
				m_particleMass);
	}

	void PICPlusPlus::initializeSpecies(double& particleCharge, 
		double& particleMass,
		double& chargeCloudWidth, 
		std::vector<double>& particlePositions, 
		std::vector<double>& particleXVelocities,
		const double plasmaFrequency, 
		const int	 numParticles, 
		const double chargeMassRatio,  
		const double speciesPlasmaFrequency, 
		const double driftVelocity, 
		const double thermalVelocity) {
		
		particleCharge = m_spatialLength * plasmaFrequency * plasmaFrequency / (numParticles * chargeMassRatio);
		particleMass = particleCharge / chargeMassRatio;
		chargeCloudWidth = m_spatialLength / numParticles;

		initializeLinearPositions(particlePositions, numParticles, chargeCloudWidth);

		addInitialVelocities(particleXVelocities, numParticles, driftVelocity, thermalVelocity);

		for (int K = 0; K < numParticles; ++K) {
			particleXVelocities[K] *= m_dtdx;
		}
	}

	void PICPlusPlus::initializeLinearPositions(std::vector<double>& inOutParticlePositions, const int numParticles, const double chargeCloudWidth) {
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

	void PICPlusPlus::addInitialVelocities(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity, const double thermalVelocity) {
		if (driftVelocity != 0) {
			addDriftVelocity(inOutParticleXVelocities, numParticles, driftVelocity);
		}

		if (thermalVelocity != 0) {
			addThermalVelocity(inOutParticleXVelocities, numParticles, thermalVelocity);
		}
	}

	void PICPlusPlus::applySpatialPerturbation(std::vector<double>& inOutParticlePositions, const int numParticles, const int spatialPerturbationMode, const double spatialPerturbationAmplitude) {
		for (int a = 0; a < numParticles; ++a) {
			double theta = 2 * std::numbers::pi * spatialPerturbationMode * inOutParticlePositions[a] / m_spatialLength;
			inOutParticlePositions[a] = inOutParticlePositions[a] + spatialPerturbationAmplitude * cos(theta);

			if (inOutParticlePositions[a] >= m_spatialLength) {
				inOutParticlePositions[a] = inOutParticlePositions[a] - m_spatialLength;
			}
			if (inOutParticlePositions[a] < 0) {
				inOutParticlePositions[a] = inOutParticlePositions[a] + m_spatialLength;
			}
		}
	}

	void PICPlusPlus::calculateEnergies(
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass) {

		for (int species = 0; species < m_numSpecies; species++) {
			for (int i = 0; i < m_speciesNumParticles[species]; i++) {
				m_particleKineticEnergy[species][m_timeStep] += 0.5 * std::pow((particleXVelocities[species][i]), 2) * particleMass[species];
			}
		}

		for (int i = 0; i < m_numGrid + 1; i++) {
			m_electrostaticEnergy[m_timeStep] += std::pow(m_electricField[m_timeStep][i], 2) * 0.5 * m_gridStepSize;
		}
	}

	std::vector<DATA_STRUCTS::Particle> PICPlusPlus::updateFrameParticles(
		int numSpecies,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particlePositions,
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass
	) {
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