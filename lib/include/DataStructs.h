#pragma once

#include <nlohmann/json.hpp>

namespace DATA_STRUCTS {

	struct SimulationParams {
		int numGrid;
		double spatialLength;
		double gridStepSize;
		int numTimeSteps;
		double timeStepSize;
		int numSpecies;
	};

	struct SpeciesData {
		std::string name;
		int numParticles;
		double driftVelocity;
		double thermalVelocity;
		double spatialPerturbationAmplitude;
		int spatialPerturbationMode;
		int plasmaFrequency;
		int chargeMassRatio;
		double particleCharge;
		double particleMass;
		double chargeCloudWidth;

		std::vector<double> particlePositions;
		std::vector<double> particleXVelocities;
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
