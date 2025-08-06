cmake_minimum_required(VERSION 3.15)
project(emerald VERSION 2.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify binary output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Create library target
add_library(emerald_lib
    src/core/config.cpp
    src/io/cli_handler.cpp
    src/analysis/alignment_processor.cpp
    src/io/alignment_writer.cpp
)

target_include_directories(emerald_lib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Create main executable
add_executable(emerald src/main.cpp)
target_link_libraries(emerald PRIVATE emerald_lib)

# Testing configuration
enable_testing()
find_package(GTest REQUIRED)

# Unit tests
add_executable(unit_tests
    tests/unit/core/config_test.cpp
    tests/unit/io/cli_handler_test.cpp
)

target_link_libraries(unit_tests
    PRIVATE
    emerald_lib
    GTest::GTest
    GTest::Main
)

add_test(NAME unit_tests COMMAND unit_tests)

# Create test directories
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/data/input)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/data/output)
