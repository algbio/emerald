#include "dp_matrix.h"
#include <algorithm>

DPMatrix::DPMatrix(size_t n, size_t m) : n_(n), m_(m) {
    initializeMatrices();
}

void DPMatrix::initializeMatrices() {
    edge_matrix_.assign(n_ + 1, 
        std::vector<std::vector<std::vector<AlignmentEdge>>>(m_ + 1,
        std::vector<std::vector<AlignmentEdge>>(3)));
}

void DPMatrix::buildAdjacencyMatrix(const std::string& seq_a, const std::string& seq_b, 
                                   const ScoringMatrix& scoring, const Config& config) {
    for (size_t i = 0; i <= n_; ++i) {
        for (size_t j = 0; j <= m_; ++j) {
            addTransitions(i, j, seq_a, seq_b, scoring, config);
        }
    }
}

void DPMatrix::addTransitions(size_t i, size_t j, const std::string& seq_a, 
                             const std::string& seq_b, const ScoringMatrix& scoring, 
                             const Config& config) {
    // Gap in sequence B (horizontal move) - edge from (i,j) to (i,j+1)
    if (j + 1 <= m_) {
        edge_matrix_[i][j][stateToIndex(AlignmentState::MATCH)].emplace_back(
            i, j + 1, AlignmentState::GAP_B, config.start_gap + config.gap_cost);
        edge_matrix_[i][j][stateToIndex(AlignmentState::GAP_B)].emplace_back(
            i, j + 1, AlignmentState::GAP_B, config.gap_cost);
    }
    
    // Gap in sequence A (vertical move) - edge from (i,j) to (i+1,j)
    if (i + 1 <= n_) {
        edge_matrix_[i][j][stateToIndex(AlignmentState::MATCH)].emplace_back(
            i + 1, j, AlignmentState::GAP_A, config.start_gap + config.gap_cost);
        edge_matrix_[i][j][stateToIndex(AlignmentState::GAP_A)].emplace_back(
            i + 1, j, AlignmentState::GAP_A, config.gap_cost);
    }
    
    // Match/mismatch (diagonal move) - edge from (i,j) to (i+1,j+1)
    if (i + 1 <= n_ && j + 1 <= m_) {
        int64_t match_cost = scoring.getScore(seq_a[i], seq_b[j]);
        edge_matrix_[i][j][stateToIndex(AlignmentState::MATCH)].emplace_back(
            i + 1, j + 1, AlignmentState::MATCH, match_cost);
    }
    
    // State transitions (no position move, just state change)
    edge_matrix_[i][j][stateToIndex(AlignmentState::GAP_A)].emplace_back(
        i, j, AlignmentState::MATCH, 0);
    edge_matrix_[i][j][stateToIndex(AlignmentState::GAP_B)].emplace_back(
        i, j, AlignmentState::MATCH, 0);
}

DPMatrix::ScoreMatrix DPMatrix::computeForwardDP() {
    ScoreMatrix dp(n_ + 1, 
        std::vector<std::vector<int64_t>>(m_ + 1,
        std::vector<int64_t>(3, -(1LL << 30))));
    
    dp[0][0][stateToIndex(AlignmentState::MATCH)] = 0;
    
    for (size_t i = 0; i <= n_; ++i) {
        for (size_t j = 0; j <= m_; ++j) {
            for (int k = 2; k >= 0; --k) {
                if (dp[i][j][k] <= -(1LL << 30)) continue;
                
                for (const auto& edge : edge_matrix_[i][j][k]) {
                    int target_state = stateToIndex(edge.getTargetState());
                    dp[edge.getTargetN()][edge.getTargetM()][target_state] = std::max(
                        dp[edge.getTargetN()][edge.getTargetM()][target_state], 
                        dp[i][j][k] + edge.getTransitionCost()
                    );
                }
            }
        }
    }
    
    return dp;
}

DPMatrix::ScoreMatrix DPMatrix::computeBackwardDP() {
    ScoreMatrix dpr(n_ + 1, 
        std::vector<std::vector<int64_t>>(m_ + 1,
        std::vector<int64_t>(3, -(1LL << 30))));
    
    dpr[n_][m_][stateToIndex(AlignmentState::MATCH)] = 0;
    
    for (int64_t i = n_; i >= 0; --i) {
        for (int64_t j = m_; j >= 0; --j) {
            for (int k = 0; k <= 2; ++k) {
                for (const auto& edge : edge_matrix_[i][j][k]) {
                    int target_state = stateToIndex(edge.getTargetState());
                    if (dpr[edge.getTargetN()][edge.getTargetM()][target_state] <= -(1LL << 30)) continue;
                    dpr[i][j][k] = std::max(
                        dpr[i][j][k],
                        dpr[edge.getTargetN()][edge.getTargetM()][target_state] + edge.getTransitionCost()
                    );
                }
            }
        }
    }
    
    return dpr;
}