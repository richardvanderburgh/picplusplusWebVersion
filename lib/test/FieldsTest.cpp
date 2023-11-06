#include <gtest/gtest.h>

#include <Fields.hpp>

TEST(FieldsTest, fields) {
	std::vector<double> rho(33, 1.0);
	double L = 1.0;
	double dx = 1.0;
	std::vector<std::vector<double>> E(2, std::vector<double>(33, 1.0));
	int t = 1;
	const int constNg = 1;
	std::vector<double> a(33, 0.0);
	double ael = 1.0;
	fields(rho, L, dx, E, t, constNg, a, ael);

	//EXPECT_EQ(rho, std::vector<double>());
	//EXPECT_EQ(E, std::vector<std::vector<double>>());
	//EXPECT_EQ(a, std::vector<double>());
	EXPECT_EQ(ael, 1.0);
}
