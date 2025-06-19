#ifndef WASM_INTERFACE_H
#define WASM_INTERFACE_H

#include <string>
#include <cstdint>

// C++ style function declaration
std::string generate_alignment_json(const std::string& representative_sequence,
                                  const std::string& representative_descriptor,
                                  const std::string& member_sequence,
                                  const std::string& member_descriptor,
                                  float alpha = 0.75f,
                                  int64_t delta = 0,
                                  int64_t gap_cost = -1,
                                  int64_t start_gap = -11);

// C-style function for backward compatibility
extern "C" {
#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#define C_EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define C_EXPORT
#endif

C_EXPORT const char* generate_alignment_json_c(const char* representative_sequence,
                                             const char* representative_descriptor,
                                             const char* member_sequence,
                                             const char* member_descriptor,
                                             float alpha,
                                             int64_t delta,
                                             int64_t gap_cost,
                                             int64_t start_gap);
}

#endif // WASM_INTERFACE_H