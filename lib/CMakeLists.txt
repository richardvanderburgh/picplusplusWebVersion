cmake_minimum_required(VERSION 3.5)
project(PIC++Lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# find_package(Matplot++ REQUIRED HINTS $ENV{MATPLOTPP_DIR})

add_library(${PROJECT_NAME} STATIC

include/Accel.hpp
include/First_EE.hpp
include/Fields.hpp
include/init.hpp
include/fft.hpp
include/complex.hpp
include/SetRho.hpp
include/Utils.hpp

fft.cpp
complex.cpp
)

#include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
#conan_basic_setup()


target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries( ${PROJECT_NAME} 
    #${CONAN_LIBS}
    nlohmann_json::nlohmann_json
    # Matplot++::matplot
)

