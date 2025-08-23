#include "alignment_processor.h"
#include "dag_builder.h"
#include "safety_window_calculator.h"
#include <iostream>

void AlignmentProcessor::processSequences(const std::vector<std::unique_ptr<Sequence>>& sequences) {
    if (config_.reference_protein.empty()) {
        processAllPairs(sequences);
    } else {
        processWithReference(sequences);
    }
}

void AlignmentProcessor::processWithReference(const std::vector<std::unique_ptr<Sequence>>& sequences) {
    size_t ref_idx = findReferenceSequence(sequences);
    
    for (size_t i = 0; i < sequences.size(); ++i) {
        if (i != ref_idx) {
            computeAlignment(*sequences[ref_idx], *sequences[i]);
        }
    }
}

void AlignmentProcessor::processAllPairs(const std::vector<std::unique_ptr<Sequence>>& sequences) {
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (size_t j = i + 1; j < sequences.size(); ++j) {
            computeAlignment(*sequences[i], *sequences[j]);
        }
    }
}

size_t AlignmentProcessor::findReferenceSequence(const std::vector<std::unique_ptr<Sequence>>& sequences) {
    for (size_t i = 0; i < sequences.size(); ++i) {
        if (sequences[i]->getDescriptor() == config_.reference_protein) {
            return i;
        }
    }
    throw std::runtime_error("Reference protein not found: " + config_.reference_protein);
}

void AlignmentProcessor::computeAlignment(const Sequence& seq1, const Sequence& seq2) {
    if (seq1.getAlphabet()->getType() != seq2.getAlphabet()->getType()) {
        throw std::runtime_error("Cannot align sequences of different types: " + 
                               seq1.getAlphabet()->getName() + " vs " + 
                               seq2.getAlphabet()->getName());
    }
    
    // Ensure scoring matrix matches sequence type
    if (scoring_matrix_.getAlphabetType() != seq1.getAlphabet()->getType()) {
        throw std::runtime_error("Scoring matrix type doesn't match sequence type");
    }
    
    const std::string& a = seq1.getSequence();
    const std::string& b = seq2.getSequence();
    
    validateSequences(a, seq1.getDescriptor());
    validateSequences(b, seq2.getDescriptor());
    
    if (config_.verbose) {
        std::cout << "Aligning: " << seq1.getDescriptor() 
                  << " vs " << seq2.getDescriptor() 
                  << " (lengths: " << a.size() << " x " << b.size() << ")" << std::endl;
    }
    
    auto result = computeOptimalAlignment(a, b);
    alignment_writer_.writeAlignment(seq1, seq2, result);
}

AlignmentProcessor::AlignmentResult 
AlignmentProcessor::computeOptimalAlignment(const std::string& a, const std::string& b) {
    const size_t n = a.size();
    const size_t m = b.size();
    
    // Create DP matrix and compute alignment
    DPMatrix dp_matrix(n, m);
    dp_matrix.buildAdjacencyMatrix(a, b, scoring_matrix_, config_);
    
    auto forward_scores = dp_matrix.computeForwardDP();
    auto backward_scores = dp_matrix.computeBackwardDP();
    
    const int64_t optimal_score = forward_scores[n][m][0];  // Changed to int64_t
    
    if (config_.verbose) {
        std::cout << "Optimal alignment score: " << optimal_score << std::endl;
    }
    
    // Build delta-neighborhood DAG
    DAGBuilder dag_builder(config_);
    auto alignment_result = dag_builder.buildDAG(dp_matrix.getEdgeMatrix(),
                               forward_scores, backward_scores, 
                               optimal_score, n, m);

    // Compute safety windows using the DAG and DP scores
    SafetyWindowCalculator sw_calc(config_);
    auto windows = sw_calc.computeSafetyWindows(
        alignment_result, forward_scores, backward_scores, n, m);

    // Attach to result (expects AlignmentResult to have `std::vector<SafetyWindow> safety_windows`)
    alignment_result.safety_windows = std::move(windows);

    return alignment_result;
}

void AlignmentProcessor::validateSequences(const std::string& seq, const std::string& name) const {
    for (size_t i = 0; i < seq.size(); ++i) {
        if (!scoring_matrix_.isValidCharacter(seq[i])) {
            throw std::runtime_error("Invalid character '" + std::string(1, seq[i]) + 
                                   "' at position " + std::to_string(i) + " in " + name);
        }
    }
}