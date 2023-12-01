#pragma once

#include <optional>
#include <vector>

#include <DataStructs.h>

namespace PIC_PLUS_PLUS {

	class PICPlusPlus {

	public:
		PICPlusPlus(DATA_STRUCTS::InputVariables inputVariables);


		[[nodiscard]] std::optional<nlohmann::json> initialize();

	private:
		DATA_STRUCTS::SimulationParams m_simulationParams;
		std::vector<DATA_STRUCTS::SpeciesData> m_allSpeciesData;

		std::vector<double> m_chargeCloudWidth;
		std::vector<double> m_qdx;

		std::vector<std::vector<double>> m_rhos;
		std::vector<std::vector<double>> m_rho0;

		double m_dtdx;
		double m_ael;

		int m_timeStep;

		std::vector<double> m_particleAcceleration;
		std::vector<double> m_chargeDensity;
		std::vector<double> m_electrostaticEnergy;
		std::vector<double> m_totalEnergy;

		std::vector<std::vector<double>> m_electricField;
		std::vector<std::vector<double>> m_particleKineticEnergy;
		std::vector<std::vector<double>> m_particleDriftEnergy;

		struct PicData {
			std::vector<DATA_STRUCTS::Frame> frames;
		};

		PicData mPicData;

		void runTimeLoop();

		void initializePositions(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const double chargeCloudWidth);

		void addDriftVelocity(std::vector<double>& inOutParticleXVelocities,
			const int numParticles, 
			const double driftVelocity);

		void addThermalVelocity(std::vector<double>& inOutParticleXVelocities, 
			const int numParticles, 
			const double thermalVelocity);

		void initializeVelocities(std::vector<double>& inOutParticleXVelocities, 
			const int numParticles, 
			const double driftVelocity, 
			const double thermalVelocity);

		void applySpatialPerturbation(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const int spatialPerturbationMode,
			const double spatialPerturbationAmplitude);

		void calculateEnergies();

		DATA_STRUCTS::Frame updateFrame();

		[[nodiscard]] std::vector<DATA_STRUCTS::Particle> updateFrameParticles();
	};
}
