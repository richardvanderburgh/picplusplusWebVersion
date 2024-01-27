add_executable(${PROJECT_NAME}_Test
    LoggerTest.cpp
)

target_include_directories(${PROJECT_NAME}_Test
    PRIVATE
    ./../private_headers
)

target_link_libraries(${PROJECT_NAME}_Test
    PRIVATE
    ${PROJECT_NAME}
    GTest::gtest_main
)

add_test(NAME "${PROJECT_NAME}"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
    COMMAND ${PROJECT_NAME}_Test)