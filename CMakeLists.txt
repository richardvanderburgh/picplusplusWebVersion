cmake_minimum_required(VERSION 3.5)
project(PIC++Main LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

enable_testing()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

find_package(GTest REQUIRED)
#find_package(FFTW3 REQUIRED)
find_package(nlohmann_json REQUIRED)

add_subdirectory(lib)

add_executable(${PROJECT_NAME}
    src/main.cpp
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
    PIC++Lib
)

add_executable(${PROJECT_NAME}_Test
    test/PICplusplusTest.cpp
    test/AccelTest.cpp
    test/FieldsTest.cpp
    test/PICTest.cpp
    test/UtilsTest.cpp
)

target_link_libraries(${PROJECT_NAME}_Test
    PRIVATE
    PIC++Lib
    GTest::gtest
)

add_test(NAME "${PROJECT_NAME}"
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
          COMMAND ${PROJECT_NAME}_Test)
