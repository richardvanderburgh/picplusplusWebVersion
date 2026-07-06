cmake_minimum_required(VERSION 3.5)
project(PIC++Lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(PUBLIC_HEADERS
	./include/DataStructs.h
	./include/PICPlusPlus.h
)

set(PRIVATE_HEADERS
    ./private_headers/Accel.hpp
    ./private_headers/Fields.hpp
    ./private_headers/fft.hpp
    ./private_headers/complex.hpp
    ./private_headers/SetRho.hpp
    ./private_headers/Utils.hpp
)

set(SOURCE_FILES
	./src/complex.cpp
	./src/fft.cpp
    ./src/PICPlusPlus.cpp
)

add_library(${PROJECT_NAME} STATIC
	${PUBLIC_HEADERS}
	${PRIVATE_HEADERS}
	${SOURCE_FILES}
)

target_include_directories(${PROJECT_NAME}
	PUBLIC
	  ${CMAKE_CURRENT_SOURCE_DIR}/include
	PRIVATE
	  ${CMAKE_CURRENT_SOURCE_DIR}/private_headers
)

target_link_libraries( ${PROJECT_NAME} 
    PUBLIC
    nlohmann_json::nlohmann_json
    PRIVATE
    project_warnings
)

# Propagate OpenMP (compile flags + runtime) to everything that links PIC++Lib:
# the main executable and both test binaries. PUBLIC ensures the _OPENMP define
# reaches the header-only hot loops (Accel.hpp / SetRho.hpp) wherever they are
# compiled.
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

add_subdirectory(test)

