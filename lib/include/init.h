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
			int chargeMassRatio);

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
#endif // !INIT_HPP