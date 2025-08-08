// This file is a stripped-down version for implementing the extern "C" functions
// without Emscripten-specific includes that cause compilation issues

#include "wasm_interface.h"
#include <string>

// C-style wrappers implemented separately to avoid header conflicts
extern "C" {

const char* generate_alignment_json_c(const char* refSeq, const char* refDesc, 
                                    const char* memSeq, const char* memDesc,
                                    float alpha, int64_t delta, int64_t gapCost, int64_t startGap) {
    static std::string result_str;
    
    // Call the C++ version
    result_str = generate_alignment_json(refSeq, refDesc, memSeq, memDesc, 
                                      alpha, delta, gapCost, startGap);
    
    // Return pointer to the static buffer (will persist until next call)
    return result_str.c_str();
}

const char* generate_alignment_json_with_matrix_c(const char* refSeq, const char* refDesc,
                                               const char* memSeq, const char* memDesc,
                                               int matrix_type, float alpha, int64_t delta, 
                                               int64_t gapCost, int64_t startGap) {
    static std::string result_str;
    
    // Call the C++ version with the selected matrix type
    result_str = generate_alignment_json_with_matrix(refSeq, refDesc, memSeq, memDesc,
                                                  static_cast<CostMatrixType>(matrix_type),
                                                  alpha, delta, gapCost, startGap);
    
    // Return pointer to the static buffer (will persist until next call)
    return result_str.c_str();
}

} // end extern "C"
