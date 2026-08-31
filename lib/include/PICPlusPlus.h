#pragma once

#include <optional>
#include <vector>

#include <DataStructs.h>
#include "Grid3D.hpp"

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

		bool m_is3D = false;
		Grid3D m_grid3D;
		std::vector<double> m_chargeDensity3D;
		std::vector<double> m_electricFieldX3D;
		std::vector<double> m_electricFieldY3D;
		std::vector<double> m_electricFieldZ3D;

		struct PicData {
			std::vector<DATA_STRUCTS::Frame> frames;
		};

		PicData mPicData;

		void runTimeLoop();
		void runTimeLoop3D();

		void initializePositions(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const double chargeCloudWidth);

		void initializePositions3D(DATA_STRUCTS::SpeciesData& speciesData);

		void addDriftVelocity(std::vector<double>& inOutParticleXVelocities,
			const int numParticles, 
			const double driftVelocity);

		void addThermalVelocity(std::vector<double>& inOutParticleXVelocities, 
			const int numParticles, 
			const double thermalVelocity);

		void addThermalVelocity3D(DATA_STRUCTS::SpeciesData& speciesData);

		void initializeVelocities(std::vector<double>& inOutParticleXVelocities, 
			const int numParticles, 
			const double driftVelocity, 
			const double thermalVelocity);

		void initializeVelocities3D(DATA_STRUCTS::SpeciesData& speciesData);

		void applySpatialPerturbation(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const int spatialPerturbationMode,
			const double spatialPerturbationAmplitude,
			const std::string& spatialPerturbationWaveform);

		void applySpatialPerturbation3D(DATA_STRUCTS::SpeciesData& speciesData);

		void calculateEnergies();
		void calculateEnergies3D();

		DATA_STRUCTS::Frame updateFrame();
		DATA_STRUCTS::Frame updateFrame3D();

		[[nodiscard]] std::vector<DATA_STRUCTS::Particle> updateFrameParticles();
		[[nodiscard]] std::vector<DATA_STRUCTS::Particle> updateFrameParticles3D();

		[[nodiscard]] std::optional<nlohmann::json> initialize1D();
		[[nodiscard]] std::optional<nlohmann::json> initialize3D();
	};
}
