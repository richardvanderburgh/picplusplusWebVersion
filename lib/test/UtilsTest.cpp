#include <gtest/gtest.h>

#include <Utils.hpp>

TEST(UtilsTest, linespace) {
    const double start = 1.0;
    const double end = 1000000.0;
    const int num_points = 1000000;
    //const std::vector<double> truth{5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0};
    std::vector<double> truth(num_points);
    double step = (end - start) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        truth[i] = start + i * step;
    }
    
    const auto result = linspace(start, end, num_points);
    ASSERT_EQ(result, truth);
}