#include "alignment_processor.h"
#include <iostream>

void AlignmentProcessor::computeAlignment(const Sequence& seq1, const Sequence& seq2) {
    // Stub implementation that creates dummy alignment results
    std::vector<std::pair<size_t, size_t>> dummy_windows;
    dummy_windows.emplace_back(0, 10);
    dummy_windows.emplace_back(20, 30);
    
    if (config_.verbose) {
        std::cout << "Processing alignment between: \n"
                  << seq1.getDescriptor() << " and "
                  << seq2.getDescriptor() << std::endl;
    }
    
    // Write different output formats based on config
    alignment_writer_.writeAlignmentResults(dummy_windows, seq1, seq2);
    
    if (!config_.json_file.empty()) {
        alignment_writer_.writeJsonReport(dummy_windows, seq1, seq2);
    }
    
    if (config_.draw_graph) {
        alignment_writer_.writeAlignmentGraph(dummy_windows, seq1, seq2);
    }
}