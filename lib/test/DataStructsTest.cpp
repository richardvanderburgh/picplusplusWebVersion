#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <DataStructs.h>

TEST(DataStructsTest, particle_toJson) {
	DATA_STRUCTS::Particle particle{.position=2.0, .velocity=5.0, .species=2, .id=10};
	nlohmann::json jsonTruth = nlohmann::json{ {"position", particle.position}, {"velocity", particle.velocity}, {"species", particle.species}, {"id", particle.id} };

	nlohmann::json jsonParticle = particle;
	EXPECT_EQ(jsonParticle, jsonTruth);
}

TEST(DataStructsTest, particle_fromJson) {
	DATA_STRUCTS::Particle particleTruth{.position=2.0, .velocity=5.0, .species=2, .id=10};
	nlohmann::json particleJson = nlohmann::json{ {"position", particleTruth.position}, {"velocity", particleTruth.velocity}, {"species", particleTruth.species}, {"id", particleTruth.id} };

	DATA_STRUCTS::Particle particle = particleJson.get<DATA_STRUCTS::Particle>();
	EXPECT_EQ(particle, particleTruth);
}

TEST(DataStructsTest, frame_toJson) {
	DATA_STRUCTS::Frame frame{ .particles = {DATA_STRUCTS::Particle{}, DATA_STRUCTS::Particle{}}, .electricField = {1.0, 2.0}, .frameNumber = 2 };
	nlohmann::json jsonTruth = nlohmann::json{ {"particles", frame.particles}, {"electricField", frame.electricField}, {"frameNumber", frame.frameNumber} };

	nlohmann::json jsonFrame = frame;
	EXPECT_EQ(jsonFrame, jsonTruth);
}

TEST(DataStructsTest, frame_fromJson) {
	DATA_STRUCTS::Frame frameTruth{ .particles = {DATA_STRUCTS::Particle{}, DATA_STRUCTS::Particle{}}, .electricField = {1.0, 2.0}, .frameNumber = 2 };
	nlohmann::json frameJson = nlohmann::json{ {"particles", frameTruth.particles}, {"electricField", frameTruth.electricField}, {"frameNumber", frameTruth.frameNumber} };

	DATA_STRUCTS::Frame frame = frameJson.get<DATA_STRUCTS::Frame>();
	EXPECT_EQ(frame, frameTruth);
}