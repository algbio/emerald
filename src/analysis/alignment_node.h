#pragma once
#include "alignment_state.h"
#include <cstdint>
#include <functional>

class AlignmentNode {
public:
    using State = AlignmentState;  // Use shared enum
    
    AlignmentNode(size_t i, size_t j, State state)
        : i_(static_cast<size_t>(i)), 
          j_(static_cast<size_t>(j)), 
          state_(state) {}
    
    // Default constructor for containers
    AlignmentNode() : i_(0), j_(0), state_(State::MATCH) {}
    
    size_t getI() const { return i_; }
    size_t getJ() const { return j_; }
    State getState() const { return state_; }
    
    // Conversion to unique ID (for compatibility with existing code)
    size_t toId(size_t m) const { return i_ * (m + 1) + j_; }
    
    // Factory method from ID
    static AlignmentNode fromId(size_t id, size_t m, State state) {
        return AlignmentNode(id / (m + 1), id % (m + 1), state);
    }
    
    // Comparison operators for use in containers
    bool operator==(const AlignmentNode& other) const {
        return i_ == other.i_ && j_ == other.j_ && state_ == other.state_;
    }
    
    bool operator<(const AlignmentNode& other) const {
        if (i_ != other.i_) return i_ < other.i_;
        if (j_ != other.j_) return j_ < other.j_;
        return state_ < other.state_;
    }
    
private:
    size_t i_, j_;
    State state_;
};
