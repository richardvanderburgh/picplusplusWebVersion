project(Logger LANGUAGES CXX)

set(PUBLIC_HEADERS
    ./include/Logger.h
)

set(PRIVATE_HEADERS
)

set(SOURCE_FILES
    ./src/Logger.cpp
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
    spdlog::spdlog
)

add_subdirectory(test)