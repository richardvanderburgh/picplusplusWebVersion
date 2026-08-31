#pragma once

#include <nlohmann/json.hpp>

namespace DATA_STRUCTS {

	struct SimulationParams {
		int dimension = 1;
		int numGrid;
		int numGridY = 0;
		int numGridZ = 0;
		double spatialLength;
		double spatialLengthY = 0.0;
		double spatialLengthZ = 0.0;
		double gridStepSize;
		double gridStepSizeY = 0.0;
		double gridStepSizeZ = 0.0;
		int numTimeSteps;
		double timeStepSize;
		int numSpecies;
		// Record a phase-space/field frame every `framePeriod` steps. 1 (default)
		// records every step; values <= 0 disable per-step frame recording, which
		// is useful for large performance runs where storing every frame would
		// dominate runtime and memory.
		int framePeriod = 1;
		// Uniform external magnetic field (electrostatic PIC + Lorentz force).
		// When any component is nonzero, the Boris pusher is used and 1D runs
		// become 1D3V (three velocity components on a 1D mesh).
		double magneticFieldX = 0.0;
		double magneticFieldY = 0.0;
		double magneticFieldZ = 0.0;

		[[nodiscard]] bool hasMagneticField() const {
			return magneticFieldX != 0.0 || magneticFieldY != 0.0 || magneticFieldZ != 0.0;
		}

		[[nodiscard]] bool usesVelocity3V() const {
			return dimension == 3 || hasMagneticField();
		}
	};

	struct SpeciesData {
		std::string name;
		int numParticles;
		double driftVelocity;
		double driftVelocityY = 0.0;
		double driftVelocityZ = 0.0;
		double thermalVelocity;
		double spatialPerturbationAmplitude;
		int spatialPerturbationMode;
		// "cos" (default) or "sin". A sin shift seeds a cosine density mode and is the
		// standard standing-wave initial condition at t = 0.
		std::string spatialPerturbationWaveform = "cos";
		double plasmaFrequency;
		double chargeMassRatio;
		double particleCharge;
		double particleMass;
		double chargeCloudWidth;

		std::vector<double> particlePositions;
		std::vector<double> particlePositionsY;
		std::vector<double> particlePositionsZ;
		std::vector<double> particleXVelocities;
		std::vector<double> particleYVelocities;
		std::vector<double> particleZVelocities;
	};

	struct InputVariables {
		SimulationParams simulationParams;
		std::vector<SpeciesData> allSpeciesData;
	};


	struct Particle {
		double position = 0.0;
		double positionY = 0.0;
		double positionZ = 0.0;
		double velocity = 0.0;
		double velocityY = 0.0;
		double velocityZ = 0.0;
		int32_t species = 0;
		int32_t id = 0;

		[[nodiscard]] bool operator==(const Particle&) const noexcept = default;
	};

	NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Particle, position, positionY, positionZ, velocity, velocityY, velocityZ, species, id);

	struct Frame {
		std::vector<Particle> particles{};
		std::vector<double> electricField{};
		int32_t frameNumber = 0;
		int32_t dimension = 1;

		[[nodiscard]] bool operator==(const Frame&) const noexcept = default;
	};

	NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(Frame, particles, electricField, frameNumber, dimension);
}
