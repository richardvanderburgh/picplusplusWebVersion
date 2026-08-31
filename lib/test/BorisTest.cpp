#include <gtest/gtest.h>

#include <cmath>
#include <Boris.hpp>

TEST(BorisTest, RotatePreservesSpeedAndMatchesCyclotron) {
	double vx = 0.0;
	double vy = 1.0;
	double vz = 0.0;
	const double speed0 = std::hypot(vx, vy, vz);

	// Boris angle φ satisfies tan(φ/2) = |t|. For a 90° rotation about z, |t| = 1.
	borisRotate(vx, vy, vz, 0.0, 0.0, 1.0);

	EXPECT_NEAR(std::hypot(vx, vy, vz), speed0, 1e-12);
	EXPECT_NEAR(vx, 1.0, 1e-12);
	EXPECT_NEAR(vy, 0.0, 1e-12);
	EXPECT_NEAR(vz, 0.0, 1e-12);
}

TEST(BorisTest, ZeroBReducesToElectricKick) {
	double vx = 0.0;
	double vy = 0.0;
	double vz = 0.0;
	borisPushCodeUnits(
		vx, vy, vz,
		2.0, 0.0, 0.0,
		0.5, 0.5, 0.5,
		-1.0, 0.1,
		1.0, 1.0, 1.0,
		0.0, 0.0, 0.0);
	// Total electric Δv_code = ae * Ex = 0.5 * 2 = 1
	EXPECT_NEAR(vx, 1.0, 1e-12);
	EXPECT_NEAR(vy, 0.0, 1e-12);
	EXPECT_NEAR(vz, 0.0, 1e-12);
}
