#include <unordered_map>
#include <vector>
#include <gmpxx.h>

#include "optimal_paths.h"

// Save the ratios of the path in a vector for the safety_windows function
std::vector<mpq_class> find_ratios(std::vector<int64_t> &path, std::vector<std::vector<int64_t>> &dag,
		std::vector<std::vector<mpq_class>> &ratios);

// Calculate the safety windows
std::tuple<std::vector<std::pair<int64_t, int64_t>>, std::vector<mpq_class>, int64_t>
safety_windows(std::vector<mpz_class> &am, std::vector<mpz_class> &ram,
		std::vector<int64_t> &path, mpq_class alpha);

// Prints out edges that are safe but not optimal, motivating the exploration
// of the suboptimal space in a sequence alignment (only run in verbose mode)
std::vector<std::pair<int64_t, int64_t>> safe_not_opt(std::vector<int64_t> &path,
		        std::vector<std::pair<int64_t, int64_t>> &swindows, Dag &d, bool for_representative);
