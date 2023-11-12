#pragma once

#include <nlohmann/json.hpp>

namespace DATA_STRUCTS {

	struct SimulationConfig {
		double spatialLength;
		int numTimeSteps;
		double timeStepSize;
		int numGrid;
	};

	struct SpeciesProperties {
		double plasmaFrequency = 0;
		double chargeMassRatio = 0;
		double thermalVelocity = 0;
		double driftVelocity = 0;
		int numParticles = 0;
		int spatialPerturbationMode = 0;
		double spatialPerturbationAmplitude = 0;
	};

	struct Particle {
		double position = 0.0;
		double velocity = 0.0;
		int32_t species = 0;
		int32_t id = 0;

		[[nodiscard]] bool operator==(const Particle&) const noexcept = default;
	};

	NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Particle, position, velocity, species, id);

	struct Frame {
		std::vector<Particle> particles{};
		std::vector<double> electricField{};
		int32_t frameNumber = 0;

		[[nodiscard]] bool operator==(const Frame&) const noexcept = default;
	};

	NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Frame, particles, electricField, frameNumber);
}
