//#include <init.hpp>
#include <gtest/gtest.h>
#include "include/PICTest.hpp"


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
    return 0;
}

