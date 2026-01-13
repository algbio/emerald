#include "dag_builder.h"
#include "alignment_node.h"
#include "graph_edge.h"
#include "alignment_result.h"
#include "safety_window_calculator.h"
#include <iostream>

DAGBuilder::AlignmentResult DAGBuilder::buildDAG(
    const DPMatrix::EdgeMatrix& edge_matrix,
    const DPMatrix::ScoreMatrix& forward_scores,
    const DPMatrix::ScoreMatrix& backward_scores,
    int64_t optimal_score, size_t n, size_t m) {
    
    AlignmentResult result(optimal_score);
    int64_t delta_threshold = optimal_score - config_.delta;
    
    // Pre-allocate space for efficiency
    result.delta_neighborhood.reserve(3*n + 3*m);
    
    if (config_.verbose) {
        std::cout << "Building delta-neighborhood graph (optimal score: " << optimal_score 
                  << ", delta: " << config_.delta 
                  << ", threshold: " << delta_threshold << ")\n";
    }
    
    // Build delta-neighborhood DAG - all alignments within delta of optimal
    for (size_t i = 0; i <= n; ++i) {
        for (size_t j = 0; j <= m; ++j) {
            for (int k = 0; k <= 2; ++k) {
                if (forward_scores[i][j][k] <= -(1LL << 30)) continue;
                
                // Create source node once for this state
                AlignmentNode source_node(i, j, static_cast<AlignmentState>(k));
                
                for (const auto& edge : edge_matrix[i][j][k]) {
                    int target_state_idx = static_cast<int>(edge.getTargetState());
                    int64_t path_score = forward_scores[i][j][k] + edge.getTransitionCost() + 
                                       backward_scores[edge.getTargetN()][edge.getTargetM()][target_state_idx];
                    
                    //std::cerr << path_score << ' ' << delta_threshold << std::endl;
                    //std::cerr << forward_scores[i][j][k] << ' ' << backward_scores[edge.getTargetN()][edge.getTargetM()][target_state_idx] << ' ' << edge.getTransitionCost() << std::endl;

                    // Include edge in delta-neighborhood if within threshold
                    if (path_score >= delta_threshold) {
                        AlignmentNode target_node(
                            edge.getTargetN(), 
                            edge.getTargetM(), 
                            edge.getTargetState()
                        );
                        result.delta_neighborhood.emplace_back(
                            source_node, 
                            target_node, 
                            edge.getTransitionCost()
                        );
                        std::cerr << "Adding edge: " << result.delta_neighborhood.back().toString() << " with path score " << path_score << std::endl;
                    }
                }
            }
        }
    }
    
    if (config_.verbose) {
        std::cout << "Delta-neighborhood contains " << result.delta_neighborhood.size() << " edges\n";
    }
    
    // Calculate safety windows from the delta neighborhood
    SafetyWindowCalculator window_calculator(config_);
    result.safety_windows = window_calculator.computeSafetyWindows(
        result, forward_scores, backward_scores, n, m);
    
    // Compute statistics
    result.statistics["neighborhood_size"] = static_cast<double>(result.delta_neighborhood.size());
    result.statistics["delta_threshold"] = static_cast<double>(delta_threshold);
    result.statistics["sequence_lengths"] = static_cast<double>(n * m);
    result.statistics["coverage_ratio"] = static_cast<double>(result.delta_neighborhood.size()) / (n * m);
    result.statistics["num_safety_windows"] = static_cast<double>(result.safety_windows.size());
    
    if (config_.verbose) {
        std::cout << "Found " << result.safety_windows.size() << " safety windows\n";
    }
    
    return result;
}