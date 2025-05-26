#include <vector>
#include <string>
#include <gmpxx.h>

#include "optimal_paths.h"

void write_json_file(const std::string &json_file, Protein &ref, Protein &mem, Dag &d, std::vector<std::pair<int64_t, int64_t>> &windows, std::vector<std::pair<int64_t, int64_t>> &windowsp, std::vector<std::vector<mpq_class>> &ratio);
