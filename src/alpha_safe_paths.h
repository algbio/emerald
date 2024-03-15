#pragma once
#include <vector>
#include <gmpxx.h>

#include "optimal_paths.h"

// Given a dag of optimal paths, find a path with almost safe (>= alpha) paths
std::vector<int64_t> topsort(std::vector<std::vector<int64_t>> &dag);

// Find (any) path from src to dest. Assertion error in case if such path does not exist
void find_path(int64_t src, int64_t dest, std::vector<int64_t> &path,
		std::vector<std::vector<int64_t>> &dag, std::vector<int64_t> &order);

// For each vertex, save the number of paths starting from the vertex to the sink
std::vector<mpz_class> number_of_paths(std::vector<std::vector<int64_t>> &dag);

// For each edge, calculate the percentage of s-t paths they are part in
std::vector<std::vector<mpq_class>> path_ratios(Dag &d, std::vector<mpz_class> &am, std::vector<mpz_class> &ram);

// Find an s to t path that contains all edges with occurence ratio >= alpha. Might fail if alpha <= 0.5
std::vector<int64_t> find_alpha_path(Dag &d,
		std::vector<std::vector<mpq_class>> &ratios, mpq_class alpha, bool verbose_flag);

