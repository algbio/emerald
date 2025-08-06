#include "config.h"
#include <filesystem>
#include <iostream>

void Config::validateConfig() const {
    validateFiles();
    validateParameters();
}

void Config::validateFiles() const {
    namespace fs = std::filesystem;
    
    if (input_file.empty()) {
        throw std::runtime_error("Input file must be specified");
    }
    if (!fs::exists(input_file)) {
        throw std::runtime_error("Input file does not exist: " + input_file);
    }
    
    if (output_file.empty()) {
        throw std::runtime_error("Output file must be specified");
    }
    
    if (!cost_matrix_file.empty() && !fs::exists(cost_matrix_file)) {
        throw std::runtime_error("Cost matrix file does not exist: " + cost_matrix_file);
    }
    
    if (draw_graph && !fs::exists(draw_graph_dir)) {
        throw std::runtime_error("Graph output directory does not exist: " + draw_graph_dir);
    }
}

void Config::validateParameters() const {
    if (alpha <= 0.5f || alpha > 1.0f) {
        throw std::runtime_error("Alpha must be in range (0.5, 1.0]");
    }
    
    if (gap_cost >= 0) {
        throw std::runtime_error("Gap cost must be negative");
    }
    
    if (start_gap >= 0) {
        throw std::runtime_error("Start gap cost must be negative");
    }
}