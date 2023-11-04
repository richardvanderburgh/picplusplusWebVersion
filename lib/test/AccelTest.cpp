#include <gtest/gtest.h>

#include <Accel.hpp>

TEST(AccelTest, accel) {
	const int nsp = 5;
	double dx = 1.0;
	double dt = 1.0;
	int t = 0;
	std::vector<double> q(nsp, 1.0);
	std::vector<double> m(nsp, 1.0);
	double ael = 1.0;
	int ng = 1;
	std::vector<double> a(ng+1, 1.0);
	std::vector<int> N(nsp, 1.0);
	std::vector <std::vector<double>> x(nsp, std::vector<double>(1.0));
	std::vector <std::vector<double>> vx(nsp, std::vector<double>(1.0));
	accel(nsp, dx, dt, t, q, m, ael, a, ng, N, x, vx);

	EXPECT_EQ(ael, -0.5);
	EXPECT_EQ(a, std::vector<double>({-0.5, -0.5}));
	EXPECT_EQ(x, std::vector<std::vector<double>>({ { 0 }, { 0 }, { 0 }, { 0 }, { 0 } }));
	EXPECT_EQ(vx, std::vector<std::vector<double>>({ { -0.5 }, { -0.5 }, { -0.5 }, { -0.5 }, { -0.5 } }));
}