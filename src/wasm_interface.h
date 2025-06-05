#ifndef WASM_INTERFACE_H
#define WASM_INTERFACE_H

#include <string>
#include <vector>
#include <gmpxx.h>

// We'll use the existing structures rather than creating new ones
// Just forward declare them here if needed
struct Protein;
struct Dag;

#ifdef __cplusplus
extern "C" {
#endif

// Main WASM-exposed function that generates JSON from input sequences
// We need extern "C" to prevent name mangling
const char* generate_alignment_json(const char* representative_sequence,
                                  const char* representative_descriptor,
                                  const char* member_sequence,
                                  const char* member_descriptor,
                                  float alpha,
                                  int64_t delta,
                                  int64_t gap_cost,
                                  int64_t start_gap);

#ifdef __cplusplus
}
#endif

#endif // WASM_INTERFACE_H