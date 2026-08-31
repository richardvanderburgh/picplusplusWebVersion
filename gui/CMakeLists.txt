find_package(Qt6 REQUIRED COMPONENTS Widgets Charts DataVisualization)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

add_executable(PIC++GUI
    main.cpp
    MainWindow.cpp
    ChartPanel.cpp
    ParticleView3D.cpp
    SimulationConfig.cpp
    SimulationWorker.cpp
)

target_include_directories(PIC++GUI
    PRIVATE
    ${CMAKE_SOURCE_DIR}/lib/private_headers
    ${CMAKE_SOURCE_DIR}/lib/include
)

target_link_libraries(PIC++GUI
    PRIVATE
    PIC++Lib
    Qt6::Widgets
    Qt6::Charts
    Qt6::DataVisualization
    project_warnings
)

set_target_properties(PIC++GUI PROPERTIES
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
)
