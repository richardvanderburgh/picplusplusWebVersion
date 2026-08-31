cmake_minimum_required(VERSION 3.5)
project(PIC++Main LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

enable_testing()

option(PICPP_ENABLE_WARNINGS "Enable a strict compiler warning set for project code" ON)
option(PICPP_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)
option(PICPP_ENABLE_OPENMP "Enable OpenMP shared-memory parallelism when available" ON)
option(PICPP_BUILD_GUI "Build the Qt desktop GUI (requires Qt 6 Widgets, Charts, and DataVisualization)" ON)

# OpenMP is optional: when found, the hot particle loops run multi-threaded and
# the compiler defines _OPENMP (which the pragmas in the headers key off of).
# When not found (e.g. stock AppleClang without libomp), the pragmas compile out
# and the code runs correctly in serial.
if(PICPP_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX)

    # AppleClang does not ship OpenMP; Homebrew's keg-only libomp needs hints.
    if(NOT OpenMP_CXX_FOUND AND APPLE)
        foreach(_libomp_root IN ITEMS "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
            if(EXISTS "${_libomp_root}/include/omp.h")
                message(STATUS "PIC++: probing Homebrew libomp at ${_libomp_root}")
                set(OpenMP_ROOT "${_libomp_root}")
                set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I${_libomp_root}/include")
                set(OpenMP_CXX_LIB_NAMES omp)
                set(OpenMP_omp_LIBRARY "${_libomp_root}/lib/libomp.dylib")
                find_package(OpenMP COMPONENTS CXX)
                if(OpenMP_CXX_FOUND)
                    break()
                endif()
            endif()
        endforeach()
    endif()

    if(OpenMP_CXX_FOUND)
        message(STATUS "PIC++: OpenMP enabled (${OpenMP_CXX_VERSION})")
    else()
        message(STATUS "PIC++: OpenMP not found - building serial")
    endif()
endif()

# INTERFACE target carrying the warning flags. Linked only into project code
# (not third-party dependencies like GTest / nlohmann_json), so their warnings
# don't drown out ours.
add_library(project_warnings INTERFACE)
if(PICPP_ENABLE_WARNINGS)
    target_compile_options(project_warnings INTERFACE
        $<$<CXX_COMPILER_ID:MSVC>:/W4 $<$<BOOL:${PICPP_WARNINGS_AS_ERRORS}>:/WX>>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-Wall;-Wextra;-Wpedantic;-Wshadow;-Wconversion;-Wsign-conversion;-Wnon-virtual-dtor;-Wold-style-cast;-Wcast-align;-Wunused;-Woverloaded-virtual;-Wnull-dereference;-Wdouble-promotion;$<$<BOOL:${PICPP_WARNINGS_AS_ERRORS}>:-Werror>>
    )
endif()

find_package(GTest REQUIRED)
find_package(nlohmann_json REQUIRED)

add_subdirectory(lib)

add_executable(${PROJECT_NAME}
    src/main.cpp
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
    lib/private_headers
)

target_link_libraries(${PROJECT_NAME}
    PUBLIC
    PIC++Lib
    PRIVATE
    project_warnings
)

add_executable(${PROJECT_NAME}_Test
    test/PICplusplusTest.cpp
    test/PICTest.cpp
    test/ValidationTest.cpp
)

target_include_directories(${PROJECT_NAME}_Test
    PRIVATE
    lib/private_headers
)

target_link_libraries(${PROJECT_NAME}_Test
    PRIVATE
    PIC++Lib
    GTest::gtest
    project_warnings
)

add_test(NAME "${PROJECT_NAME}"
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
          COMMAND ${PROJECT_NAME}_Test)

if(PICPP_BUILD_GUI)
    # Homebrew Qt is keg-only; help CMake find it on macOS when CMAKE_PREFIX_PATH
    # was not set by the user.
    if(APPLE AND NOT DEFINED Qt6_DIR)
        foreach(_qt_root IN ITEMS "/opt/homebrew/opt/qt" "/usr/local/opt/qt"
                                  "/opt/homebrew/opt/qt@6" "/usr/local/opt/qt@6")
            if(EXISTS "${_qt_root}/lib/cmake/Qt6")
                list(PREPEND CMAKE_PREFIX_PATH "${_qt_root}")
                message(STATUS "PIC++: probing Homebrew Qt at ${_qt_root}")
                break()
            endif()
        endforeach()
    elseif(UNIX AND NOT APPLE AND NOT DEFINED Qt6_DIR)
        foreach(_qt_root IN ITEMS "/usr/lib/x86_64-linux-gnu/cmake" "/usr/lib/cmake")
            if(EXISTS "${_qt_root}/Qt6/Qt6Config.cmake")
                list(PREPEND CMAKE_PREFIX_PATH "${_qt_root}")
                message(STATUS "PIC++: probing system Qt at ${_qt_root}")
                break()
            endif()
        endforeach()
    endif()

    find_package(Qt6 QUIET COMPONENTS Widgets Charts DataVisualization)
    if(Qt6Widgets_FOUND AND Qt6Charts_FOUND AND Qt6DataVisualization_FOUND)
        add_subdirectory(gui)
        message(STATUS "PIC++: Qt GUI enabled (Qt ${Qt6_VERSION})")
    else()
        message(STATUS "PIC++: Qt6 Widgets/Charts/DataVisualization not found - skipping PIC++GUI")
        if(APPLE)
            message(STATUS "       Install with: brew install qt")
            message(STATUS "       Or set: export CMAKE_PREFIX_PATH=\"\$(brew --prefix qt)\"")
        else()
            message(STATUS "       On Debian/Ubuntu: sudo apt install qt6-base-dev qt6-charts-dev qt6-datavisualization-dev qt6-declarative-dev")
        endif()
    endif()
endif()
