#include "alignment_writer.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>

AlignmentWriter::AlignmentWriter(const Config& config) : config_(config) {}

void AlignmentWriter::writeAlignment(
    const Sequence& seq1, 
    const Sequence& seq2, 
    const AlignmentResult& result) {
    
    std::ofstream outfile;
    if (!config_.output_file.empty()) {
        outfile.open(config_.output_file);
        if (!outfile.is_open()) {
            throw std::runtime_error("Cannot open output file: " + config_.output_file);
        }
    }
    
    // Use either file or std::cout
    std::ostream& out = config_.output_file.empty() ? std::cout : outfile;
    
    out << "# Alignment Results\n";
    out << "- Sequence 1: " << seq1.getDescriptor() << " (" << seq1.getSequence().size() << " residues)\n";
    out << "- Sequence 2: " << seq2.getDescriptor() << " (" << seq2.getSequence().size() << " residues)\n";
    out << "- Optimal score: " << result.optimal_score << "\n";
    out << "- Delta neighborhood size: " << result.delta_neighborhood.size() << " edges\n\n";
    
    // Display safety windows
    out << "# Safety Windows\n";
    for (size_t i = 0; i < result.safety_windows.size(); ++i) {
        const auto& window = result.safety_windows[i];
        out << "Window " << (i+1) << ": " 
            << "Positions " << window.start_pos << "-" << window.end_pos 
            << " (target: " << window.target_start << "-" << window.target_end << ")"
            << ", Confidence: " << std::fixed << std::setprecision(2) << window.confidence << "\n";
            
        // Extract and show the actual alignment in this window
        std::string seq1_window = seq1.getSequence().substr(window.start_pos, 
                                                           window.end_pos - window.start_pos);
        std::string seq2_window = seq2.getSequence().substr(window.target_start, 
                                                           window.target_end - window.target_start);
        
        // Show the alignment
        out << seq1_window << "\n";
        
        // Generate match indicators
        std::string match_indicators;
        size_t min_len = std::min(seq1_window.length(), seq2_window.length());
        for (size_t j = 0; j < min_len; ++j) {
            match_indicators += (seq1_window[j] == seq2_window[j]) ? '|' : ' ';
        }
        out << match_indicators << "\n";
        
        out << seq2_window << "\n\n";
    }
    
    // Show statistics
    out << "# Statistics\n";
    for (const auto& [key, value] : result.statistics) {
        out << "- " << key << ": " << value << "\n";
    }
}