#ifndef INIT_HPP
#define INIT_HPP

#include <optional>
#include <vector>

#include "DataStructs.h"

namespace PIC_PLUS_PLUS {
	class Init {
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

		void initializeLinearPositions(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const double chargeCloudWidth);

		void addDriftVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity);

		void addThermalVelocities(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double thermalVelocity);

		void applySpatialPerturbation(std::vector<double>& inOutParticlePositions, 
			const int numParticles, 
			const int spatialPerturbationMode, 
			const double spatialLength, 
			const double spatialPerturbationAmplitude);

		[[nodiscard]] std::vector<DATA_STRUCTS::Particle> updateFrameParticles(
			int numSpecies,
			const std::vector<int>& speciesNumParticles,
			const std::vector<std::vector<double>>& particlePositions,
			const std::vector<std::vector<double>>& particleXVelocities,
			std::vector<std::vector<double>>& particleKineticEnergy,
			const int timeStep,
			const std::vector<double>& particleMass
		);
	};
}
#endif