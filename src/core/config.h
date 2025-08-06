#pragma once
#include <string>
#include <cstdint>

class Config {
public:
    // Alignment parameters
    float alpha = 0.75f;
    std::uint32_t delta = 0;

    // Scoring parameters
    int gap_cost = -1;
    int start_gap = -11;
    int special_score = -1;

    // Control flags
    bool ignore_special = false;
    bool verbose = false;
    bool window_merge = false;
    bool draw_graph = false;
    
    // File paths
    std::string input_file;
    std::string output_file;
    std::string cost_matrix_file;
    std::string json_file;
    std::string reference_protein;
    std::string draw_graph_dir = ".";
    std::string alignments_file;

    // Validation methods
    void validateConfig() const;

private:
    void validateFiles() const;
    void validateParameters() const;

    friend class CliHandler;  // Allow CliHandler to modify config
};