#ifndef WASM_INTERFACE_H
#define WASM_INTERFACE_H

#include <string>
#include <cstdint>
#include <vector>
#include <utility> // for std::pair

// Include optimal_paths.h for Dag and Protein structures
#include "optimal_paths.h"

// Emscripten-specific definitions
#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define C_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define C_EXPORT
#endif

// Enum for predefined cost matrices
enum CostMatrixType {
    BLOSUM62,
    PAM250,
    IDENTITY
};

// C++ style function declarations with default cost matrix (BLOSUM62)
std::string generate_alignment_json(const std::string& representative_sequence,
                                  const std::string& representative_descriptor,
                                  const std::string& member_sequence,
                                  const std::string& member_descriptor,
                                  float alpha = 0.75f,
                                  int64_t delta = 0,
                                  int64_t gap_cost = -1,
                                  int64_t start_gap = -11);

// C++ style function with predefined matrix selection
std::string generate_alignment_json_with_matrix(const std::string& representative_sequence,
                                              const std::string& representative_descriptor,
                                              const std::string& member_sequence,
                                              const std::string& member_descriptor,
                                              CostMatrixType matrix_type,
                                              float alpha = 0.75f,
                                              int64_t delta = 0,
                                              int64_t gap_cost = -1,
                                              int64_t start_gap = -11);

// C++ style function with custom cost matrix
std::string generate_alignment_json_with_custom_matrix(const std::string& representative_sequence,
                                                     const std::string& representative_descriptor,
                                                     const std::string& member_sequence,
                                                     const std::string& member_descriptor,
                                                     const std::vector<std::vector<int64_t>>& custom_matrix,
                                                     float alpha = 0.75f,
                                                     int64_t delta = 0,
                                                     int64_t gap_cost = -1,
                                                     int64_t start_gap = -11);

// Helper function to get a predefined cost matrix
const int64_t (*get_cost_matrix(CostMatrixType matrix_type))[21];

// Helper function to generate alignment strings from a path
std::pair<std::string, std::string> generate_alignment_strings(Dag &d, const std::vector<int64_t> &path, const Protein &a, const Protein &b);

// C-style functions for backward compatibility
extern "C" {
C_EXPORT const char* generate_alignment_json_c(const char* representative_sequence,
                                             const char* representative_descriptor,
                                             const char* member_sequence,
                                             const char* member_descriptor,
                                             float alpha,
                                             int64_t delta,
                                             int64_t gap_cost,
                                             int64_t start_gap);

C_EXPORT const char* generate_alignment_json_with_matrix_c(const char* representative_sequence,
                                                         const char* representative_descriptor,
                                                         const char* member_sequence,
                                                         const char* member_descriptor,
                                                         int matrix_type,
                                                         float alpha,
                                                         int64_t delta,
                                                         int64_t gap_cost,
                                                         int64_t start_gap);
}

// Standalone function wrappers for EMSCRIPTEN_BINDINGS
std::string generateAlignmentJson(const std::string& refSeq, 
                                const std::string& refDesc,
                                const std::string& memSeq, 
                                const std::string& memDesc,
                                float alpha = 0.75, 
                                int64_t delta = 0,
                                int64_t gapCost = -1, 
                                int64_t startGap = -11);

std::string generateAlignmentJsonWithMatrix(const std::string& refSeq, 
                                          const std::string& refDesc,
                                          const std::string& memSeq, 
                                          const std::string& memDesc,
                                          int matrix_type,
                                          float alpha = 0.75, 
                                          int64_t delta = 0,
                                          int64_t gapCost = -1, 
                                          int64_t startGap = -11);

std::string generateAlignmentJsonWithCustomMatrix(const std::string& refSeq, 
                                                const std::string& refDesc,
                                                const std::string& memSeq, 
                                                const std::string& memDesc,
                                                const std::vector<std::vector<int64_t>>& custom_matrix,
                                                float alpha = 0.75, 
                                                int64_t delta = 0,
                                                int64_t gapCost = -1, 
                                                int64_t startGap = -11);

#endif // WASM_INTERFACE_H