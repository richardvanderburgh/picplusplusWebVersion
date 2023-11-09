#ifndef INIT_HPP
#define INIT_HPP

#include <optional>
#include <vector>

#include "DataStructs.h"

namespace PIC_PLUS_PLUS {
	class PICPlusPlus {
	public:

		struct PicData {
			std::vector<DATA_STRUCTS::Frame> frames;
		};

		PicData mPicData;

		[[nodiscard]] std::optional<nlohmann::json> initialize(
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
			const int chargeMassRatio);

		void initializeSpecies(double& particleCharge, double& particleMass, double& chargeCloudWidth, std::vector<double>& particlePositions, std::vector<double>& particleXVelocities,
			const double spatialLength,
			const double plasmaFrequency,
			const int numParticles,
			const double chargeMassRatio,
			const double speciesPlasmaFrequency,
			const double driftVelocity,
			const double thermalVelocity,
			const double dtdx);

		void initializeLinearPositions(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const double chargeCloudWidth);

		void addDriftVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity);

		void addThermalVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double thermalVelocity);

		void addInitialVelocities(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity, const double thermalVelocity);

		void applySpatialPerturbation(std::vector<double>& inOutParticlePositions, 
			const int numParticles, 
			const int spatialPerturbationMode, 
			const double spatialLength, 
			const double spatialPerturbationAmplitude);
		
		void calculateEnergies(
			std::vector<std::vector<double>>& inOutParticleKineticEnergy,
			std::vector<double>& inOutElectrostaticEnergy,
			const int numSpecies,
			const int timeStep,
			const std::vector<int>& speciesNumParticles,
			const std::vector<std::vector<double>>& particleXVelocities,
			const std::vector<double>& particleMass,
			const int numGrid,
			const std::vector<std::vector<double>>& electricField,
			const double gridStepSize);

		[[nodiscard]] std::vector<DATA_STRUCTS::Particle> updateFrameParticles(
			int numSpecies,
			const std::vector<int>& speciesNumParticles,
			const std::vector<std::vector<double>>& particlePositions,
			const std::vector<std::vector<double>>& particleXVelocities,
			const int timeStep,
			const std::vector<double>& particleMass
		);
		DATA_STRUCTS::Frame updateFrame(
			const std::vector<std::vector<double>>& electricField,
			const int timeStep, const int numSpecies,
			const std::vector<int>& speciesNumParticles,
			const std::vector<std::vector<double>>& particlePositions,
			const std::vector<std::vector<double>>& particleXVelocities,
			const std::vector<double>& particleMass);
	};
}
#endif