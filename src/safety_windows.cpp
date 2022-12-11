#include <iostream>
#include <map>
#include <vector>

#include "safety_windows.h"

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
