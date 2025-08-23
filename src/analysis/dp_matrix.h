#pragma once
#include "alignment_edge.h"
#include "scoring_matrix.h"
#include "../core/config.h"
#include <vector>
#include <string>

class DPMatrix {
public:
    using EdgeMatrix = std::vector<std::vector<std::vector<std::vector<AlignmentEdge>>>>;
    using ScoreMatrix = std::vector<std::vector<std::vector<int64_t>>>;  // Changed to int64_t
    
    DPMatrix(size_t n, size_t m);
    
    void buildAdjacencyMatrix(const std::string& seq_a, const std::string& seq_b, 
                             const ScoringMatrix& scoring, const Config& config);
    
    ScoreMatrix computeForwardDP();
    ScoreMatrix computeBackwardDP();
    
    const EdgeMatrix& getEdgeMatrix() const { return edge_matrix_; }
    size_t getN() const { return n_; }
    size_t getM() const { return m_; }
    
private:
    size_t n_, m_;
    EdgeMatrix edge_matrix_;  
    
    void initializeMatrices();
    void addTransitions(size_t i, size_t j, const std::string& seq_a, 
                       const std::string& seq_b, const ScoringMatrix& scoring, 
                       const Config& config);
                       
    static int stateToIndex(AlignmentState state) {
        return static_cast<int>(state);
    }
};