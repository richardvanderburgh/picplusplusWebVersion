
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
		m_chargeMassRatio(chargeMassRatio) {

		m_gridStepSize = m_spatialLength / m_numGrid;
		m_dtdx = m_timeStepSize / m_gridStepSize;
	};

	std::optional<nlohmann::json> PICPlusPlus::initialize() {

		std::vector<int> speciesNumParticles(m_numSpecies, m_numParticles);
		std::vector<double> speciesPlasmaFrequency(m_numSpecies, m_plasmaFrequency);
		std::vector<double> speciesChargeMassRatio(m_numSpecies, m_chargeMassRatio);
		std::vector<double> speciesThermalVelocity(m_numSpecies, m_thermalVelocity);

		const std::vector<double> speciesDriftVelocity = { m_driftVelocity, - m_driftVelocity, 0, 2 *  m_driftVelocity, -2 *  m_driftVelocity };

		std::vector<int>    speciesSpatialPerturbationMode(m_numSpecies, m_spatialPerturbationMode);
		std::vector<double> speciesSpatialPerturbationAmplitude(m_numSpecies, m_spatialPerturbationAmplitude);

		std::vector<std::vector<double>> particleKineticEnergy(m_numSpecies, std::vector<double>(m_numTimeSteps + 1, 0.0));
		std::vector<std::vector<double>> particleDriftEnergy(m_numSpecies, std::vector<double>(m_numTimeSteps + 1, 0.0));
		std::vector<std::vector<double>> particleThermalEnergy(m_numSpecies, std::vector<double>(m_numTimeSteps + 1, 0.0));

		int maxN = *std::max_element(speciesNumParticles.begin(), speciesNumParticles.end());
		std::vector<std::vector<double>> particlePositions(m_numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleXVelocities(m_numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleYVelocities(m_numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleZVelocities(m_numSpecies, std::vector<double>(maxN, 0.0));

		std::vector<double> electrostaticEnergy(m_numTimeSteps + 1, 0.0);
		std::vector<double> totalEnergy(m_numTimeSteps + 1, 0.0);
		std::vector<double> chargeCloudWidth(m_numSpecies);
		std::vector<double> particleCharge(m_numSpecies, 0.0);
		std::vector<double> particleMass(m_numSpecies);
		std::vector<double> qdx(m_numSpecies);

		std::vector<std::vector<double>> electricField(m_numTimeSteps + 1, std::vector<double>(m_numGrid + 1, 0.0));
		std::vector<double> chargeDensity(m_numGrid + 1, 0.0);
		std::vector<double> particleAcceleration(m_numGrid + 1);

		std::vector<std::vector<double>> rhos(m_numSpecies, std::vector<double>(m_numGrid + 1, 0.0));
		std::vector<std::vector<double>> rho0(m_numSpecies, std::vector<double>(m_numGrid + 1, 0.0));

		for (int species = 0; species < m_numSpecies; species++) {

			initializeSpecies(particleCharge[species], 
				particleMass[species], 
				chargeCloudWidth[species], 
				particlePositions[species], 
				particleXVelocities[species],
				speciesPlasmaFrequency[species],
				speciesNumParticles[species],
				speciesChargeMassRatio[species],
				speciesPlasmaFrequency[species],
				speciesDriftVelocity[species],
				speciesThermalVelocity[species]);

			if (speciesSpatialPerturbationAmplitude[species] != 0) {
				applySpatialPerturbation(particlePositions[species], 
					speciesNumParticles[species], 
					speciesSpatialPerturbationMode[species],
					speciesSpatialPerturbationAmplitude[species]);
			}

			qdx[species] = particleCharge[species] / m_gridStepSize;
			setRho(species, m_numGrid, m_gridStepSize, speciesNumParticles, qdx, chargeDensity, particlePositions, rho0, rhos);
		}

		int timeStep = 0;
		double ael = 1;

		fields(chargeDensity, m_spatialLength, m_gridStepSize, electricField, timeStep, m_numGrid, particleAcceleration, ael);

		accel(particleAcceleration,
			particleCharge,
			m_numSpecies,
			m_gridStepSize,
			m_timeStepSize,
			timeStep,
			particleMass,
			ael, m_numGrid, speciesNumParticles, particlePositions, particleXVelocities);

		calculateEnergies(
			particleKineticEnergy,
			electrostaticEnergy,
			timeStep,
			speciesNumParticles,
			particleXVelocities,
			particleMass,
			electricField);

		mPicData.frames.emplace_back() = updateFrame(
			electricField,
			timeStep, m_numSpecies,
			speciesNumParticles,
			particlePositions,
			particleXVelocities,
			particleMass);
 
		for (int timeStep = 1; timeStep <= m_numTimeSteps; timeStep++) {

			accel(particleAcceleration,
				particleCharge,
				m_numSpecies,
				m_gridStepSize,
				m_timeStepSize,
				timeStep,
				particleMass,
				ael, m_numGrid, speciesNumParticles, particlePositions, particleXVelocities);

			move(m_numSpecies, chargeDensity, rho0, qdx, speciesNumParticles, particlePositions, particleXVelocities, m_numGrid);
			fields(chargeDensity, m_spatialLength, m_gridStepSize, electricField, timeStep, m_numGrid, particleAcceleration, ael);

			calculateEnergies(
				particleKineticEnergy,
				electrostaticEnergy,
				timeStep,
				speciesNumParticles,
				particleXVelocities,
				particleMass,
				electricField);

			mPicData.frames.emplace_back() = updateFrame(
				electricField,
				timeStep, m_numSpecies,
				speciesNumParticles,
				particlePositions,
				particleXVelocities,
				particleMass);
		}

		nlohmann::json JSON;
		JSON["ke"] = particleKineticEnergy;
		JSON["ese"] = electrostaticEnergy;
		JSON["phaseFrames"] = mPicData.frames;

		//std::cout << JSON.dump() << std::endl;
		return JSON;
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
		std::vector<std::vector<double>>& inOutParticleKineticEnergy,
		std::vector<double>& inOutElectrostaticEnergy,
		const int timeStep,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass,
		const std::vector<std::vector<double>>& electricField) {

		for (int species = 0; species < m_numSpecies; species++) {
			for (int i = 0; i < speciesNumParticles[species]; i++) {
				inOutParticleKineticEnergy[species][timeStep] += 0.5 * std::pow((particleXVelocities[species][i]), 2) * particleMass[species];
			}
		}

		for (int i = 0; i < m_numGrid + 1; i++) {
			inOutElectrostaticEnergy[timeStep] += std::pow(electricField[timeStep][i], 2) * 0.5 * m_gridStepSize;
		}
	}

	std::vector<DATA_STRUCTS::Particle> PICPlusPlus::updateFrameParticles(
		int numSpecies,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particlePositions,
		const std::vector<std::vector<double>>& particleXVelocities,
		const int timeStep,
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
		const std::vector<std::vector<double>>& electricField,
		const int timeStep, const int numSpecies,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particlePositions,
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass) {

		DATA_STRUCTS::Frame frame;

		frame.electricField = electricField[timeStep];

		frame.particles = updateFrameParticles(
			numSpecies,
			speciesNumParticles,
			particlePositions,
			particleXVelocities,
			timeStep,
			particleMass
		);

		frame.frameNumber = timeStep;

		return frame;
	}
};