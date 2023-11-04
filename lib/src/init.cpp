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
	bool Init::initialize(
		double spatialLength,
		int numParticles,
		int numTimeSteps,
		double timeStepSize,
		int numGrid,
		int spatialPerturbationMode,
		double driftVelocity,
		int numSpecies,
		double spatialPerturbationAmplitude,
		double thermalVelocity,
		int plasmaFrequency,
		int chargeMassRatio) {

		int epsi = 1; // 1 over epsilon naught(F / m) Epsilon normalized to 1
		int particleShapeId = 2; // 1 for nearest grid point, 2 for linear

		//Species Input Variables
		std::vector<int> speciesNumParticles(numSpecies, numParticles);          // Number of simulation particles
		std::vector<double> speciesPlasmaFrequency(numSpecies, plasmaFrequency); // Plasma Frequency
		std::vector<double> speciesCyclotronFrequency(numSpecies);               // Cyclotron frequency
		std::vector<double> speciesChargeMassRatio(numSpecies, chargeMassRatio); // q / m charge to mass ratio(C / kg)
		std::vector<double> speciesThermalVelocity(numSpecies, thermalVelocity); // RMS thermal velocity for random velocities
		std::vector<double> speciesDriftVelocity = {

			driftVelocity, -driftVelocity, 0, 2 * driftVelocity, -2 * driftVelocity };	  // Drift velocity

		std::vector<int>    speciesSpatialPerturbationMode(numSpecies, spatialPerturbationMode);
		std::vector<double> speciesSpatialPerturbationAmplitude(numSpecies, spatialPerturbationAmplitude);

		double gridStepSize = spatialLength / numGrid;

		std::vector<std::vector<double>> electricField(numTimeSteps + 1, std::vector<double>(numGrid + 1, 0.0));

		int maxN = *std::max_element(speciesNumParticles.begin(), speciesNumParticles.end());

		std::vector<std::vector<double>> particleKineticEnergy(numSpecies, std::vector<double>(numTimeSteps + 1, 0.0));
		std::vector<std::vector<double>> particleMomentum(numTimeSteps + 1, std::vector<double>(numSpecies, 0.0));
		std::vector<std::vector<double>> particleDriftEnergy(numTimeSteps + 1, std::vector<double>(numSpecies, 0.0));
		std::vector<std::vector<double>> particleThermalEnergy(numTimeSteps + 1, std::vector<double>(numSpecies, 0.0));
		std::vector<std::vector<double>> particlePositions(numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleXVelocities(numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleYVelocities(numSpecies, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> particleZVelocities(numSpecies, std::vector<double>(maxN, 0.0));

		std::vector<double> electrostaticEnergy(numTimeSteps + 1, 0.0);
		std::vector<double> totalEnergy(numTimeSteps + 1, 0.0);
		std::vector<double> chargeCloudWidth(numSpecies);
		std::vector<double> particleCharge(numSpecies, 0.0);
		std::vector<double> particleMass(numSpecies);

		for (int species = 0; species < numSpecies; species++)
		{
			particleCharge[species] = spatialLength * speciesPlasmaFrequency[species] * speciesPlasmaFrequency[species] / (epsi * speciesNumParticles[species] * speciesChargeMassRatio[species]);
			particleMass[species] = particleCharge[species] / speciesChargeMassRatio[species];
			chargeCloudWidth[species] = spatialLength / speciesNumParticles[species];

			for (int I = 1; I < speciesNumParticles[species] + 1; I++) {
				particlePositions[species][I - 1] = (I - 0.5) * chargeCloudWidth[species];
			}

			for (int I = 0; I < speciesNumParticles[species]; I++) {
				particleXVelocities[species][I] = speciesDriftVelocity[species];
			}

			if (speciesThermalVelocity[species] != 0) {
				for (int I = 0; I < speciesNumParticles[species]; ++I) {
					double rm = 0;
					for (int ith = 0; ith < 12; ++ith) {
						rm += (double)rand() / RAND_MAX;
					}
					rm -= 6;

					if (speciesCyclotronFrequency[species] != 0) {
						particleYVelocities[species][I] += speciesThermalVelocity[species] * rm;
					}
					particleXVelocities[species][I] += speciesThermalVelocity[species] * rm;
				}
			}

			for (int a = 0; a < speciesNumParticles[species]; ++a) {
				double theta = 2 * std::numbers::pi * speciesSpatialPerturbationMode[species] * particlePositions[species][a] / spatialLength;
				particlePositions[species][a] = particlePositions[species][a] + speciesSpatialPerturbationAmplitude[species] * cos(theta);

				if (particlePositions[species][a] >= spatialLength) {
					particlePositions[species][a] = particlePositions[species][a] - spatialLength;
				}
				if (particlePositions[species][a] < 0) {
					particlePositions[species][a] = particlePositions[species][a] + spatialLength;
				}
			}
		}

		// Convert position to computer normalization and accumulate charge density.
		std::vector<double> qdx(numSpecies);
		std::vector<double> chargeDensity(numGrid + 1, 0.0);
		std::vector<std::vector<double>> rhos(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> rho0(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> qjdata(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> qjp1data(numSpecies, std::vector<double>(numGrid + 1, 0.0));
		std::vector<std::vector<double>> drhodata(numSpecies, std::vector<double>(numGrid + 1, 0.0));

		double dtdx = timeStepSize / gridStepSize;

		for (int species = 0; species < numSpecies; species++) {
			qdx[species] = particleCharge[species] / gridStepSize;
			setRho(species, numGrid, gridStepSize, speciesNumParticles, qdx, chargeDensity, particlePositions, rho0, rhos, particleShapeId);

			for (int K = 0; K < speciesNumParticles[species]; ++K) {
				particleXVelocities[species][K] *= dtdx;
			}
		}

		int timeStep = 0;
		double ael = 1;

		// Acceleration
		std::vector<double> particleAcceleration(numGrid + 1);
		fields(chargeDensity, spatialLength, particleShapeId, gridStepSize, electricField, timeStep, numGrid, particleAcceleration, ael);
		accel(numSpecies, gridStepSize, timeStepSize, timeStep, particleCharge, particleMass, ael, particleAcceleration, numGrid, speciesNumParticles, particlePositions, particleXVelocities);

		for (int species = 0; species < numSpecies; species++) {
			for (int i = 0; i < speciesNumParticles[species]; i++) {
				particleKineticEnergy[species][timeStep] += 0.5 * std::pow((particleXVelocities[species][i]), 2) * particleMass[species];
			}
		}

		for (int i = 0; i < numGrid + 1; i++) {
			electrostaticEnergy[timeStep] += std::pow(electricField[timeStep][i], 2) * 0.5 * gridStepSize;
		}

		DATA_STRUCTS::Frame& frame0 = mPicData.frames.emplace_back();
		std::vector<DATA_STRUCTS::Particle> particles;

		for (int species = 0; species < numSpecies; species++) {
			for (int i = 0; i < speciesNumParticles[species]; i++) {
				int particleId = 0;

				DATA_STRUCTS::Particle& particle = particles.emplace_back();

				particle.id = particleId;
				particle.position = particlePositions[species][i];
				particle.velocity = particleXVelocities[species][i];
				particle.species = species;

				particleId++;
			}
		}

		frame0.particles = particles;
		frame0.electricField = electricField[timeStep];
		frame0.frameNumber = timeStep;

		//BEGIN TIME LOOP 

		nlohmann::json JSON;
		nlohmann::json frames;

		for (int timeStep = 1; timeStep <= numTimeSteps; timeStep++) {
			DATA_STRUCTS::Frame frame;

			accel(numSpecies, gridStepSize, timeStepSize, timeStep, particleCharge, particleMass, ael, particleAcceleration, numGrid, speciesNumParticles, particlePositions, particleXVelocities);
			move(numSpecies, chargeDensity, rho0, qdx, speciesNumParticles, particlePositions, particleXVelocities, numGrid);
			fields(chargeDensity, spatialLength, particleShapeId, gridStepSize, electricField, timeStep, numGrid, particleAcceleration, ael);

			for (int i = 0; i < numGrid + 1; i++) {
				electrostaticEnergy[timeStep] += std::pow(electricField[timeStep][i], 2) * 0.5 * gridStepSize;
			}

			frame.electricField = electricField[timeStep];
			std::vector<DATA_STRUCTS::Particle> particles;

			nlohmann::json JSONFrame;
			nlohmann::json JSONParticles;
			int particleId = 0;

			for (int species = 0; species < numSpecies; species++) {
				for (int i = 0; i < speciesNumParticles[species]; i++) {

					DATA_STRUCTS::Particle& particle = particles.emplace_back();

					particle.id = particleId;
					particle.position = particlePositions[species][i];
					particle.velocity = particleXVelocities[species][i];
					particleKineticEnergy[species][timeStep] += 0.5 * std::pow((particle.velocity), 2) * particleMass[species];
					particle.species = species;

					nlohmann::json particleObject = particle;

					//particleObject["position"] = particlePositions[species][i];
					//particleObject["velocity"] = particleXVelocities[species][i];
					//particleObject["species"] = species;
					//particleObject["id"] = particleId;

					JSONParticles.push_back(particleObject);

					particleId++;
				}
			}

			JSONFrame["particles"] = JSONParticles;
			JSONFrame["frameNumber"] = timeStep;
			frames.push_back(JSONFrame);

			frame.particles = particles;
			frame.frameNumber = timeStep;
			mPicData.frames.push_back(frame);
		}

		JSON["ke"] = particleKineticEnergy;
		JSON["ese"] = electrostaticEnergy;
		JSON["phaseFrames"] = frames;

		std::cout << JSON.dump() << std::endl;

		return true;
	}
}