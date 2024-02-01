cmake_minimum_required(VERSION 3.5)
project(PIC++Lib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(PUBLIC_HEADERS
	./include/DataStructs.h
	./include/PICPlusPlus.h
    ./include/PicWindow.h
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
	./src/PicWindow.cpp
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
    nlohmann_json::nlohmann_json
	SDL2::SDL2
)

add_subdirectory(test)

