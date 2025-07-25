#ifndef WRITE_JSON_H
#define WRITE_JSON_H

#include <vector>
#include <string>
#include <gmpxx.h>
#include <iostream>

// Include headers with definitions instead of forward declarations
#include "alpha_safe_paths.h"  // Should contain Dag definition
#include "optimal_paths.h"     // Should contain Protein definition

// If the above headers don't contain the definitions,
// you'll need to include whichever headers actually define these structures

void write_json_file(const std::string &json_file, Protein &ref, Protein &mem, Dag &d, 
                    std::vector<std::pair<int64_t, int64_t>> &windows, 
                    std::vector<std::pair<int64_t, int64_t>> &windowsp, 
                    std::vector<std::vector<mpq_class>> &ratio);

// Overload with alignment strings
void write_json_file(const std::string &json_file, Protein &ref, Protein &mem, Dag &d, 
                    std::vector<std::pair<int64_t, int64_t>> &windows, 
                    std::vector<std::pair<int64_t, int64_t>> &windowsp, 
                    std::vector<std::vector<mpq_class>> &ratio,
                    const std::string &alignment_ref,
                    const std::string &alignment_mem);

// Function to write to a stream instead of file
void write_json_to_stream(std::ostream &output_stream, Protein &ref, Protein &mem, Dag &d, 
                         std::vector<std::pair<int64_t, int64_t>> &windows, 
                         std::vector<std::pair<int64_t, int64_t>> &windowsp, 
                         std::vector<std::vector<mpq_class>> &ratio,
                         const std::string &alignment_ref,
                         const std::string &alignment_mem);

#endif // WRITE_JSON_H
