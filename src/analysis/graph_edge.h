#pragma once
#include "alignment_node.h"
#include <cstdint>
#include <string>

struct GraphEdge {
    AlignmentNode source;
    AlignmentNode target;
    int16_t cost;
    
    GraphEdge(const AlignmentNode& src, const AlignmentNode& tgt, int64_t c)
        : source(src), target(tgt), cost(c) {}
    
    // For debugging/output
    std::string toString() const {
        return "(" + std::to_string(source.getI()) + "," + std::to_string(source.getJ()) + "," + 
               std::to_string(static_cast<int>(source.getState())) + ") -> (" +
               std::to_string(target.getI()) + "," + std::to_string(target.getJ()) + "," + 
               std::to_string(static_cast<int>(target.getState())) + ") [cost=" + std::to_string(cost) + "]";
    }
};