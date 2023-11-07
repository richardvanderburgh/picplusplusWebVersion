#include "init.h"

#define _USE_MATH_DEFINES

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
	std::optional<nlohmann::json> Init::initialize(
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
		const int chargeMassRatio) {

		double gridStepSize = spatialLength / numGrid;
		double dtdx = timeStepSize / gridStepSize;

		std::vector<int> speciesNumParticles(numSpecies, numParticles);
		std::vector<double> speciesPlasmaFrequency(numSpecies, plasmaFrequency);
		std::vector<double> speciesChargeMassRatio(numSpecies, chargeMassRatio);
		std::vector<double> speciesThermalVelocity(numSpecies, thermalVelocity);

		const std::vector<double> speciesDriftVelocity = {
			driftVelocity, -driftVelocity, 0, 2 * driftVelocity, -2 * driftVelocity };

		std::vector<int>    speciesSpatialPerturbationMode(numSpecies, spatialPerturbationMode);
		std::vector<double> speciesSpatialPerturbationAmplitude(numSpecies, spatialPerturbationAmplitude);

		std::vector<std::vector<double>> particleKineticEnergy(numSpecies, std::vector<double>(numTimeSteps + 1, 0.0));
		std::vector<std::vector<double>> particleDriftEnergy(numSpecies, std::vector<double>(numTimeSteps + 1, 0.0));
		std::vector<std::vector<double>> particleThermalEnergy(numSpecies, std::vector<double>(numTimeSteps + 1, 0.0));

		int maxN = *std::max_element(speciesNumParticles.begin(), speciesNumParticles.end());
		std::vector<std::vector<double>> particlePositions(numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleXVelocities(numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleYVelocities(numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleZVelocities(numSpecies, std::vector<double>(maxN, 0.0));

		std::vector<double> electrostaticEnergy(numTimeSteps + 1, 0.0);
		std::vector<double> totalEnergy(numTimeSteps + 1, 0.0);
		std::vector<double> chargeCloudWidth(numSpecies);
		std::vector<double> particleCharge(numSpecies, 0.0);
		std::vector<double> particleMass(numSpecies);
		std::vector<double> qdx(numSpecies);

		std::vector<std::vector<double>> electricField(numTimeSteps + 1, std::vector<double>(numGrid + 1, 0.0));
		std::vector<double> chargeDensity(numGrid + 1, 0.0);
		std::vector<double> particleAcceleration(numGrid + 1);

		std::vector<std::vector<double>> rhos(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> rho0(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> qjdata(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> qjp1data(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> drhodata(numSpecies, std::vector<double>(numGrid + 1, 0.0));

		for (int species = 0; species < numSpecies; species++) {

			initializeSpecies(particleCharge[species], particleMass[species], chargeCloudWidth[species], particlePositions[species], particleXVelocities[species],
				spatialLength,
				speciesPlasmaFrequency[species],
				speciesNumParticles[species],
				speciesChargeMassRatio[species],
				speciesPlasmaFrequency[species],
				speciesDriftVelocity[species],
				speciesThermalVelocity[species],
				dtdx);

			if (speciesSpatialPerturbationAmplitude[species] != 0) {
				applySpatialPerturbation(particlePositions[species], speciesNumParticles[species], speciesSpatialPerturbationMode[species], spatialLength, speciesSpatialPerturbationAmplitude[species]);
			}

			qdx[species] = particleCharge[species] / gridStepSize;
			setRho(species, numGrid, gridStepSize, speciesNumParticles, qdx, chargeDensity, particlePositions, rho0, rhos);
		}

		int timeStep = 0;
		double ael = 1;

		fields(chargeDensity, spatialLength, gridStepSize, electricField, timeStep, numGrid, particleAcceleration, ael);

		accel(particleAcceleration,
			particleCharge,
			numSpecies,
			gridStepSize,
			timeStepSize,
			timeStep,
			particleMass,
			ael, numGrid, speciesNumParticles, particlePositions, particleXVelocities);

		calculateEnergies(
			particleKineticEnergy,
			electrostaticEnergy,
			numSpecies,
			timeStep,
			speciesNumParticles,
			particleXVelocities,
			particleMass,
			numGrid,
			electricField,
			gridStepSize);

		mPicData.frames.emplace_back() = updateFrame(
			electricField,
			timeStep, numSpecies,
			speciesNumParticles,
			particlePositions,
			particleXVelocities,
			particleMass);
 
		for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {

			accel(particleAcceleration,
				particleCharge,
				numSpecies,
				gridStepSize,
				timeStepSize,
				timeStep,
				particleMass,
				ael, numGrid, speciesNumParticles, particlePositions, particleXVelocities);

			move(numSpecies, chargeDensity, rho0, qdx, speciesNumParticles, particlePositions, particleXVelocities, numGrid);
			fields(chargeDensity, spatialLength, gridStepSize, electricField, timeStep, numGrid, particleAcceleration, ael);

			calculateEnergies(
				particleKineticEnergy,
				electrostaticEnergy,
				numSpecies,
				timeStep,
				speciesNumParticles,
				particleXVelocities,
				particleMass,
				numGrid,
				electricField,
				gridStepSize);

			mPicData.frames.emplace_back() = updateFrame(
				electricField,
				timeStep, numSpecies,
				speciesNumParticles,
				particlePositions,
				particleXVelocities,
				particleMass);
		}

		nlohmann::json JSON;
		JSON["ke"] = particleKineticEnergy;
		JSON["ese"] = electrostaticEnergy;
		JSON["phaseFrames"] = mPicData.frames;

		std::cout << JSON.dump() << std::endl;
		return JSON;
	}

	void Init::initializeSpecies(double& particleCharge, double& particleMass, double& chargeCloudWidth, std::vector<double>& particlePositions, std::vector<double>& particleXVelocities,
		const double spatialLength, 
		const double plasmaFrequency, 
		const int numParticles, 
		const double chargeMassRatio,  
		const double speciesPlasmaFrequency, 
		const double driftVelocity, 
		const double thermalVelocity, 
		const double dtdx) {
		
		particleCharge = spatialLength * plasmaFrequency * plasmaFrequency / (numParticles * chargeMassRatio);
		particleMass = particleCharge / chargeMassRatio;
		chargeCloudWidth = spatialLength / numParticles;

		initializeLinearPositions(particlePositions, numParticles, chargeCloudWidth);

		addInitialVelocities(particleXVelocities, numParticles, driftVelocity, thermalVelocity);

		for (int K = 0; K < numParticles; ++K) {
			particleXVelocities[K] *= dtdx;
		}
	}

	void Init::initializeLinearPositions(std::vector<double>& inOutParticlePositions, const int numParticles, const double chargeCloudWidth) {
		for (int I = 1; I < numParticles + 1; I++) {
			inOutParticlePositions[I - 1] = (I - 0.5) * chargeCloudWidth;
		}
	}

	void Init::addDriftVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity) {
		for (int I = 0; I < numParticles; I++) {
			inOutParticleXVelocities[I] = driftVelocity;
		}
	}

	void Init::addThermalVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double thermalVelocity) {

		for (int I = 0; I < numParticles; ++I) {
			double rm = 0;
			for (int ith = 0; ith < 12; ++ith) {
				rm += (double)rand() / RAND_MAX;
			}
			rm -= 6;

			inOutParticleXVelocities[I] += thermalVelocity * rm;
		}
	}

	void Init::addInitialVelocities(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity, const double thermalVelocity) {
		if (driftVelocity != 0) {
			addDriftVelocity(inOutParticleXVelocities, numParticles, driftVelocity);
		}

		if (thermalVelocity != 0) {
			addThermalVelocity(inOutParticleXVelocities, numParticles, thermalVelocity);
		}
	}

	void Init::applySpatialPerturbation(std::vector<double>& inOutParticlePositions, const int numParticles, const int spatialPerturbationMode, const double spatialLength, const double spatialPerturbationAmplitude) {
		for (int a = 0; a < numParticles; ++a) {
			double theta = 2 * std::numbers::pi * spatialPerturbationMode * inOutParticlePositions[a] / spatialLength;
			inOutParticlePositions[a] = inOutParticlePositions[a] + spatialPerturbationAmplitude * cos(theta);

			if (inOutParticlePositions[a] >= spatialLength) {
				inOutParticlePositions[a] = inOutParticlePositions[a] - spatialLength;
			}
			if (inOutParticlePositions[a] < 0) {
				inOutParticlePositions[a] = inOutParticlePositions[a] + spatialLength;
			}
		}
	}

	void Init::calculateEnergies(
		std::vector<std::vector<double>>& inOutParticleKineticEnergy,
		std::vector<double>& inOutElectrostaticEnergy,
		const int numSpecies,
		const int timeStep,
		const std::vector<int>& speciesNumParticles,
		const std::vector<std::vector<double>>& particleXVelocities,
		const std::vector<double>& particleMass,
		const int numGrid,
		const std::vector<std::vector<double>>& electricField,
		const double gridStepSize) {

		for (int species = 0; species < numSpecies; species++) {
			for (int i = 0; i < speciesNumParticles[species]; i++) {
				inOutParticleKineticEnergy[species][timeStep] += 0.5 * std::pow((particleXVelocities[species][i]), 2) * particleMass[species];
			}
		}

		for (int i = 0; i < numGrid + 1; i++) {
			inOutElectrostaticEnergy[timeStep] += std::pow(electricField[timeStep][i], 2) * 0.5 * gridStepSize;
		}

	}

	std::vector<DATA_STRUCTS::Particle> Init::updateFrameParticles(
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

	DATA_STRUCTS::Frame Init::updateFrame(
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