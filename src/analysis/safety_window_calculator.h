#pragma once
#include "alignment_result.h"
#include "dp_matrix.h"
#include "../core/config.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>

/**
 * @brief Calculates safety windows in alignments
 * 
 * A safety window represents a region of the alignment that appears
 * in a significant proportion of all near-optimal alignments.
 * These regions have high confidence and are considered reliable.
 */
class SafetyWindowCalculator {
public:
    explicit SafetyWindowCalculator(const Config& config) : config_(config) {}
    
    /**
     * @brief Computes safety windows in the alignment
     * 
     * Identifies alignment regions that appear in at least alpha fraction
     * of all source-to-sink paths in the delta neighborhood. These regions
     * are considered high-confidence parts of the alignment.
     * 
     * @param alignment_result The alignment result containing the delta neighborhood
     * @param forward_scores Forward DP scores
     * @param backward_scores Backward DP scores
     * @param n Length of first sequence
     * @param m Length of second sequence
     * @return Vector of safety windows
     */
    std::vector<SafetyWindow> computeSafetyWindows(
        const AlignmentResult& alignment_result,
        const DPMatrix::ScoreMatrix& forward_scores,
        const DPMatrix::ScoreMatrix& backward_scores,
        size_t n, size_t m);
    
    std::vector<size_t> findAlignmentWithAllSafeEdges(
        const std::vector<std::pair<size_t, size_t>> high_freq_edges,
        size_t n, size_t m);
    
private:
    const Config& config_;
};