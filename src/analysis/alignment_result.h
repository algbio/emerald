#pragma once
#include "graph_edge.h"
#include <vector>
#include <map>
#include <string>
#include <cstdint>

// Define a safety window representing high-confidence alignment regions
struct SafetyWindow {
    size_t start_pos;      // Start position in reference sequence
    size_t end_pos;        // End position in reference sequence
    double ratio;         // Confidence score (0.0-1.0)
    
    SafetyWindow(size_t start, size_t end, double ratio)
        : start_pos(start), end_pos(end),
          ratio(ratio) {}
};

struct AlignmentResult {
    int64_t optimal_score;
    std::vector<GraphEdge> delta_neighborhood;
    std::vector<SafetyWindow> safety_windows;
    std::map<std::string, double> statistics;
    
    AlignmentResult(int64_t score = 0) : optimal_score(score) {}
    
    // Helper methods
    size_t getNeighborhoodSize() const { return delta_neighborhood.size(); }
    size_t getNumWindows() const { return safety_windows.size(); }
    
    // Extract optimal path from the neighborhood if needed
    std::vector<GraphEdge> extractOptimalPath() const;
    
    // Convert to integer pairs if needed
    std::vector<std::pair<int64_t, int64_t>> getGraphAsIds(int64_t m) const {
        std::vector<std::pair<int64_t, int64_t>> ids;
        ids.reserve(delta_neighborhood.size());
        
        for (const auto& edge : delta_neighborhood) {
            ids.emplace_back(edge.source.toId(m), edge.target.toId(m));
        }
        
        return ids;
    }
};