#include <iostream>
#include <map>
#include <vector>

#include "safety_windows.h"

std::vector<mpq_class> find_ratios(std::vector<int64_t> &path, std::vector<std::vector<int64_t>> &dag,
		std::vector<std::vector<mpq_class>> &ratios) {
	int64_t k = (int64_t) path.size();
	std::vector<mpq_class> rat_path(k - 1);

	for (int64_t i = 0; i < k - 1; i++) {
		for (int64_t j = 0; j < (int64_t) dag[path[i]].size(); j++) if (dag[path[i]][j] == path[i + 1]) {
			rat_path[i] = ratios[path[i]][j];
		}
	}

	return rat_path;
}

std::tuple<std::vector<std::pair<int64_t, int64_t>>, std::vector<mpq_class>, int64_t>
safety_windows(std::vector<mpz_class> &am, std::vector<mpz_class> &ram,
		std::vector<int64_t> &path, mpq_class alpha) {
	int64_t k = (int64_t) path.size();
	if (k == 0) return {};

	int64_t n = (int64_t) am.size();

	std::unordered_map<int64_t, int64_t> order;
	for (int64_t i = 0; i < k; i++) order[path[i]] = i;

	std::vector<std::pair<int64_t, int64_t>> windows;
	std::vector<mpq_class> window_ratios;
	mpq_class a = 1;
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
		mpq_class a = 1;
		if (a * am[path[i + 1]] * ram[path[i]] / am[path[0]] >= alpha)
			number_of_edges++;
	}

	return std::make_tuple(windows, window_ratios, number_of_edges);
}

std::vector<std::pair<int64_t, int64_t>> safe_not_opt(std::vector<int64_t> &path,
		std::vector<std::pair<int64_t, int64_t>> &swindows, Dag &d, bool for_representative) {
	int64_t n = (int64_t) swindows.size();
	int64_t k = (int64_t) path.size();

	std::vector<std::pair<int64_t, int64_t>> ret;

	std::map<int64_t, int64_t> idx;
	for (int64_t i = 0; i < k; i++) idx[path[i]] = i;

	auto is_opt = [&](int64_t L, int64_t R) {
		return d.in_optimal[std::make_pair(L, R)];
	};

	int64_t current_window = 0;
	for (int64_t i = 0; i < k - 1; i++) {
		while (current_window < n && idx[swindows[current_window].second] < idx[path[i + 1]])
			current_window++;
		if (current_window == n) break;
		if (idx[path[i]] < idx[swindows[current_window].first]) continue;

		// safe edge
		assert(idx[path[i]] >= idx[swindows[current_window].first] &&
				idx[path[i + 1]] <= idx[swindows[current_window].second]);
		
		if (!is_internal(d, path[i], path[i + 1]) && !is_gap(d, path[i], path[i + 1], for_representative) && !is_opt(path[i], path[i + 1]))
			ret.emplace_back(path[i], path[i + 1]);
	}
	return ret;
}
