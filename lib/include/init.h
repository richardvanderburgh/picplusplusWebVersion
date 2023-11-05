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

		void initializeLinearPositions(const int numParticles,
			std::vector<double>& particlePositions,
			const double chargeCloudWidth);

		void addDriftVelocity(const int numParticles, std::vector<double>& particleXVelocities, const double driftVelocity);

		void addThermalVelocities(const int numParticles, std::vector<double>& particleXVelocities, const double thermalVelocity);

		void applySpatialPerturbation(const int numParticles,
			const int spatialPerturbationMode, 
			std::vector<double>& particlePositions, 
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