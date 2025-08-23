#pragma once
#include "alignment_state.h"
#include <cstdint>
#include <cstddef>

class AlignmentEdge {
public:
    using State = AlignmentState;
    
    AlignmentEdge(size_t target_n, size_t target_m, State target_state, int16_t transition_cost)
        : target_n_(static_cast<size_t>(target_n)), 
          target_m_(static_cast<size_t>(target_m)), 
          target_state_(target_state), 
          transition_cost_(static_cast<int16_t>(transition_cost)) {}
    
    // Only keep the descriptive function names
    size_t getTargetN() const { return target_n_; }
    size_t getTargetM() const { return target_m_; }
    State getTargetState() const { return target_state_; }
    int64_t getTransitionCost() const { return transition_cost_; }
    
private:
    size_t target_n_, target_m_;  // Target position after this transition
    State target_state_;            // Target state after this transition
    int64_t transition_cost_;       // Cost of this transition/edge
};