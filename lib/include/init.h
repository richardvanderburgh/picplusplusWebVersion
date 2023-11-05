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

		std::optional<nlohmann::json> initialize(
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
	};
}
#endif // !INIT_HPP