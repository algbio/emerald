cmake_minimum_required(VERSION 3.16)
project(emerald)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(GTest REQUIRED)

# Set output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directories
include_directories(src)

# Create library target with all source files
add_library(emerald_lib
    src/core/config.cpp
    src/core/sequence_alphabet.cpp
    src/io/cli_handler.cpp
    src/io/sequence_reader.cpp
    src/analysis/alignment_processor.cpp
    src/analysis/scoring_matrix.cpp
    src/analysis/dp_matrix.cpp
    src/analysis/dag_builder.cpp
    src/io/alignment_writer.cpp
    src/analysis/safety_window_calculator.cpp
)

# Main executable
add_executable(emerald
    src/main.cpp
)

target_link_libraries(emerald
    PRIVATE
    emerald_lib
)

# Unit tests
add_executable(unit_tests
    tests/unit/core/config_test.cpp
    tests/unit/io/cli_handler_test.cpp
    tests/unit/io/sequence_reader_test.cpp 
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
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
