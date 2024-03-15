#include <string>
#include <iostream>
#include <vector>

#include "optimal_paths.h"

std::string draw_subgraph(const int64_t IDX, const int64_t n, const int64_t m, const Dag &d, const std::vector<std::vector<mpq_class>> &ratios, mpq_class alpha, const std::string &a, const std::string &b);
