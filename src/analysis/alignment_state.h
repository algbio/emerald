#pragma once
#include <cstdint>

enum class AlignmentState : uint8_t { 
    MATCH = 0, 
    GAP_A = 1, 
    GAP_B = 2 
};