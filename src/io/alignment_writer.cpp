#include "alignment_writer.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <filesystem>

AlignmentWriter::AlignmentWriter(const Config& config) : config_(config) {}

void AlignmentWriter::writeAlignment(
    const Sequence& seq1, 
    const Sequence& seq2, 
    const AlignmentResult& result) {
    
    const int n = seq1.getSequence().size();
    const int m = seq2.getSequence().size();
    // Function to extract (i,j,state) from node ID
    auto idToNode = [m](size_t id) -> std::tuple<size_t, size_t, int> {
        int state = id % 3;
        size_t j = ((id - state) / 3) % (m+1);
        size_t i = (id - state - 3*j) / ((m+1) * 3);
        return {i, j, state};
    };

    std::ofstream outfile;
    if (!config_.output_file.empty()) {
        // open in append mode so we don't overwrite previous alignments
        std::filesystem::path out_path(config_.output_file);
        bool already_has_content = std::filesystem::exists(out_path) && std::filesystem::file_size(out_path) > 0;
        outfile.open(config_.output_file, std::ios::app);
        if (!outfile.is_open()) {
            throw std::runtime_error("Cannot open output file: " + config_.output_file);
        }
        // if file already had content, add a small separator for readability
        if (already_has_content) {
            outfile << "\n---\n\n";
        }
    }
    
    // Use either file or std::cout
    std::ostream& out = config_.output_file.empty() ? std::cout : outfile;
    
    out << "# Alignment Results\n";
    out << "- Sequence 1: " << seq1.getDescriptor() << " (" << seq1.getSequence().size() << " characters)\n";
    out << "- Sequence 2: " << seq2.getDescriptor() << " (" << seq2.getSequence().size() << " characters)\n";
    out << "- Optimal alignment score: " << result.optimal_score << "\n";
    out << "- Delta neighborhood size: " << result.delta_neighborhood.size() << " edges\n";
    out << "- Number of safety windows: " << result.safety_windows.size() << "\n\n";
    
    // Display safety windows
    out << "# Safety Windows\n";
    for (size_t i = 0; i < result.safety_windows.size(); ++i) {
        const auto& window = result.safety_windows[i];
        const auto& [start_rep, start_mem, _] = idToNode(window.start_pos);
        const auto& [end_rep, end_mem, __] = idToNode(window.end_pos);
        out << "Window " << (i+1) << ": "
            << "Representattive: " << start_rep << "-" << end_rep
            << ", Member: " << start_mem << "-" << end_mem
            << ", Ratio: " << std::fixed << std::setprecision(2) << window.ratio << "\n";
            
        // Extract and show the actual alignment in this window
        std::string seq1_window = seq1.getSequence().substr(start_rep, 
                                                           end_rep - start_rep);
        std::string seq2_window = seq2.getSequence().substr(start_mem, 
                                                           end_mem - start_mem);
        
        // Show the alignment
        out << seq1_window << "\n";
        out << seq2_window << "\n\n";
    }
    
    // Show statistics
    out << "# Statistics\n";
    for (const auto& [key, value] : result.statistics) {
        out << "- " << key << ": " << value << "\n";
    }
}