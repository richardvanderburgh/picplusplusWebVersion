#pragma once

#include <optional>
#include <vector>

#include "DataStructs.h"
#include <DataStructs.h>

namespace PIC_PLUS_PLUS {

	class PICPlusPlus {

	public:
		PICPlusPlus(
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


		[[nodiscard]] std::optional<nlohmann::json> initialize();

	private:
		const double m_spatialLength;
		const int m_numParticles;
		const int m_numTimeSteps;
		const double m_timeStepSize;
		const int m_numGrid;
		const int m_spatialPerturbationMode;
		const double m_driftVelocity;
		const int m_numSpecies;
		const double m_spatialPerturbationAmplitude;
		const double m_thermalVelocity;
		const int m_plasmaFrequency;
		const int m_chargeMassRatio;

		std::vector<double> m_chargeCloudWidth;
		std::vector<double> m_particleCharge;
		std::vector<double> m_particleMass;
		std::vector<double> m_qdx;

		std::vector<int>    m_speciesNumParticles;

		std::vector<std::vector<double>> m_particlePositions;
		std::vector<std::vector<double>> m_particleXVelocities;

		std::vector<std::vector<double>> m_rhos;
		std::vector<std::vector<double>> m_rho0;
		
		double m_gridStepSize;
		double m_dtdx;

		int m_timeStep;

		double m_ael;

		std::vector<double> m_particleAcceleration;
		std::vector<std::vector<double>> m_electricField;

		std::vector<double> m_chargeDensity;

		std::vector<DATA_STRUCTS::SpeciesProperties> m_speciesProperties;

		// Simulation state variables
		std::vector<std::vector<double>> m_particleKineticEnergy;
		std::vector<std::vector<double>> m_particleDriftEnergy;

		std::vector<double> m_electrostaticEnergy;
		std::vector<double> m_totalEnergy;

		struct PicData {
			std::vector<DATA_STRUCTS::Frame> frames;
		};

		PicData mPicData;


		void initializeSpecies(double& particleCharge, 
			double& particleMass, 
			double& chargeCloudWidth, 
			std::vector<double>& particlePositions, 
			std::vector<double>& particleXVelocities,
			const double plasmaFrequency,
			const int numParticles,
			const double chargeMassRatio,
			const double speciesPlasmaFrequency,
			const double driftVelocity,
			const double thermalVelocity);

		void runTimeLoop();

		void initializeLinearPositions(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const double chargeCloudWidth);

		void addDriftVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity);

		void addThermalVelocity(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double thermalVelocity);

		void addInitialVelocities(std::vector<double>& inOutParticleXVelocities, const int numParticles, const double driftVelocity, const double thermalVelocity);

		void applySpatialPerturbation(std::vector<double>& inOutParticlePositions,
			const int numParticles,
			const int spatialPerturbationMode,
			const double spatialPerturbationAmplitude);

		void calculateEnergies(
			const std::vector<std::vector<double>>& particleXVelocities,
			const std::vector<double>& particleMass);

		[[nodiscard]] std::vector<DATA_STRUCTS::Particle> updateFrameParticles(
			int numSpecies,
			const std::vector<int>& speciesNumParticles,
			const std::vector<std::vector<double>>& particlePositions,
			const std::vector<std::vector<double>>& particleXVelocities,
			const std::vector<double>& particleMass
		);
		DATA_STRUCTS::Frame updateFrame(const int numSpecies,
			const std::vector<int>& speciesNumParticles,
			const std::vector<std::vector<double>>& particlePositions,
			const std::vector<std::vector<double>>& particleXVelocities,
			const std::vector<double>& particleMass);
	};
}
