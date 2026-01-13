#pragma once
#include "sequence_alphabet.h"
#include <string>
#include <cstdint>

struct Config {
    // File paths
    std::string input_file;
    std::string output_file;
    std::string reference_protein;
    std::string cost_matrix_file;
    std::string json_file;
    
    // Sequence parameters
    SequenceType sequence_type = SequenceType::PROTEIN;
    
    // Algorithm parameters
    int64_t delta = 8;
    int64_t gap_cost = -1;
    int64_t start_gap = -11;
    int64_t special_score = -5;
    
    // Safety window parameters
    double alpha = 0.75;
    size_t min_window_size = 2;
    size_t safety_window_size = 3;
    bool window_merge = false;
    bool ignore_special = false;
    
    // Visualization options
    bool draw_graph = false;
    std::string draw_graph_dir;
    bool verbose = false;
    
    /**
     * Validate configuration parameters
     * Checks that required fields are present and parameter values are valid
     * @throws std::runtime_error if validation fails
     */
    void validateConfig() const;

    /**
     * Validate file existence
     * Checks that input files exist
     * @throws std::runtime_error if validation fails
     */
    void validateFiles() const;

    /**
     * Create necessary directories
     * Creates output directories
     * @returns true if successful, throws exception otherwise
     * @throws std::runtime_error if directory creation fails
     */
    bool createOutputDirectories();
};