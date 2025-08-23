#pragma once
#include "dp_matrix.h"
#include "alignment_result.h"
#include "../core/config.h"

class DAGBuilder {
public:
    using AlignmentResult = ::AlignmentResult;
    
    explicit DAGBuilder(const Config& config) : config_(config) {}
    
    AlignmentResult buildDAG(
        const DPMatrix::EdgeMatrix& edge_matrix,
        const DPMatrix::ScoreMatrix& forward_scores,
        const DPMatrix::ScoreMatrix& backward_scores,
        int64_t optimal_score, size_t n, size_t m);
    
private:
    const Config& config_;
};