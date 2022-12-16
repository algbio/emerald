#include <unordered_map>
#include <vector>
#include <gmpxx.h>

#include "optimal_paths.h"

// Save the ratios of the path in a vector for the safety_windows function
template<class K>
std::vector<K> find_ratios(std::vector<int64_t> &path, std::vector<std::vector<int64_t>> &dag,
		std::vector<std::vector<K>> &ratios) {
	int64_t k = (int64_t) path.size();
	std::vector<K> rat_path(k - 1);

	for (int64_t i = 0; i < k - 1; i++) {
		for (int64_t j = 0; j < (int64_t) dag[path[i]].size(); j++) if (dag[path[i]][j] == path[i + 1]) {
			rat_path[i] = ratios[path[i]][j];
		}
	}

	return rat_path;
}

// Calculate the safety windows
template<class T, class K>
std::tuple<std::vector<std::pair<int64_t, int64_t>>, std::vector<K>, int64_t> safety_windows(std::vector<T> &am, std::vector<T> &ram,
		std::vector<int64_t> &path, K alpha) {
	int64_t k = (int64_t) path.size();
	if (k == 0) return {};

	int64_t n = (int64_t) am.size();

	std::unordered_map<int64_t, int64_t> order;
	for (int64_t i = 0; i < k; i++) order[path[i]] = i;

	std::vector<std::pair<int64_t, int64_t>> windows;
	std::vector<K> window_ratios;
	K a = 1;
	auto outside = [&](const int64_t &L, const int64_t &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return order[bL] >= order[L] && order[bR] <= order[R];
	};
	auto inside = [&](const int64_t &L, const int64_t &R) {
		if (windows.empty()) return false;
		auto [bL, bR] = windows.back();
		return order[L] >= order[bL] && order[R] <= order[bR];
	};
	for (int64_t L = 0, R = 0; R < k; a = a * (R + 1 == k ? 1 : am[path[R + 1]]) / am[path[R]], R++) {
		assert((R + 1 == k ? 1 : am[path[R + 1]]) <= am[path[R]]);
		assert(path[R] < n && path[R] >= 0);
		while (L < R && a < alpha) {
			assert(ram[path[L + 1]] >= ram[path[L]]);
			a = a * ram[path[L + 1]] / ram[path[L]];
			L++;
		}
		while (outside(path[L], path[R])) windows.pop_back(), window_ratios.pop_back();
		if (L < R && !inside(path[L], path[R])) windows.emplace_back(path[L], path[R]), window_ratios.emplace_back(a);
	}
	int64_t number_of_edges = 0;
	for (int64_t i = 0; i < k - 1; i++) {
		K a = 1;
		if (a * am[path[i + 1]] * ram[path[i]] / am[path[0]] >= alpha)
			number_of_edges++;
	}

	return std::make_tuple(windows, window_ratios, number_of_edges);
}

// Prints out edges that are safe but not optimal, motivating the exploration
// of the suboptimal space in a sequence alignment (only run in verbose mode)
std::vector<std::pair<int64_t, int64_t>> safe_not_opt(std::vector<int64_t> &path,
		        std::vector<std::pair<int64_t, int64_t>> &swindows, Dag &d, bool for_representative);
