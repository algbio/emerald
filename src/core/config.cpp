#include "config.h"
#include <filesystem>
#include <iostream>

void Config::validateConfig() const {
    // Alpha parameter validation
    if (alpha <= 0.5f || alpha > 1.0f) {
        throw std::runtime_error("Alpha must be in range (0.5, 1.0]");
    }
    
    // Scoring parameters validation
    if (gap_cost >= 0) {
        throw std::runtime_error("Gap cost must be negative");
    }
    
    if (start_gap > 0) {
        throw std::runtime_error("Start gap cost can not be positive");
    }
    
    if (!ignore_special && special_score >= 0) {
        throw std::runtime_error("Special character score must be negative");
    }
    
    // Safety window parameters
    if (min_window_size < 1) {
        throw std::runtime_error("Minimum window size must be at least 1");
    }
    
    if (safety_window_size < 1) {
        throw std::runtime_error("Safety window size must be at least 1");
    }
}

void Config::validateFiles() const {
    namespace fs = std::filesystem;
    
    // Check input file existence
    if (!fs::exists(input_file)) {
        throw std::runtime_error("Input file does not exist: " + input_file);
    }
    
    // Check cost matrix file if provided
    if (!cost_matrix_file.empty() && !fs::exists(cost_matrix_file)) {
        throw std::runtime_error("Cost matrix file does not exist: " + cost_matrix_file);
    }
}

bool Config::createOutputDirectories() {
    namespace fs = std::filesystem;
    try {
        // Ensure output directory exists (if a path is provided)
        if (!output_file.empty()) {
            fs::path out_dir = fs::path(output_file).parent_path();
            if (!out_dir.empty() && !fs::exists(out_dir)) {
                fs::create_directories(out_dir);
                if (verbose) {
                    std::cout << "Created output directory: " << out_dir << std::endl;
                }
            }
        }

        if (draw_graph && !draw_graph_dir.empty()) {
            if (!fs::exists(draw_graph_dir)) {
                fs::create_directories(draw_graph_dir);
                if (verbose) {
                    std::cout << "Created graph output directory: " << draw_graph_dir << std::endl;
                }
            }
        }

        if (!json_file.empty()) {
            fs::path json_dir = fs::path(json_file).parent_path();
            if (!json_dir.empty() && !fs::exists(json_dir)) {
                fs::create_directories(json_dir);
                if (verbose) {
                    std::cout << "Created JSON output directory: " << json_dir << std::endl;
                }
            }
        }
    } catch (const fs::filesystem_error& e) {
        throw std::runtime_error("Filesystem error while creating directories: " + std::string(e.what()));
    } catch (const std::exception& e) {
        throw std::runtime_error("Error while creating directories: " + std::string(e.what()));
    }
    return true;
}