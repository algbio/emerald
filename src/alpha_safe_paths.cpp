#include <iostream>
#include <vector>
#include <set>
#include <stack>
#include <unordered_map>
#include <assert.h>

#include <gmpxx.h>

#include "alpha_safe_paths.h"
#include "optimal_paths.h"


// given a dag of optimal paths, find a path with almost safe (>= alpha) paths

std::vector<int64_t> topsort(std::vector<std::vector<int64_t>> &dag) {
	int64_t n = (int64_t) dag.size();
	std::vector<int64_t> indeg(n, 0);
	for (int64_t i = 0; i < n; i++) {
		for (int64_t v: dag[i]) indeg[v]++;
	}
	std::stack<int64_t> nxt;
	for (int64_t i = 0; i < n; i++) if (indeg[i] == 0) nxt.push(i);

	std::vector<int64_t> sorted;
	while (!nxt.empty()) {
		int64_t v = nxt.top();
		nxt.pop();
		sorted.push_back(v);

		for(int64_t u: dag[v]) if (--indeg[u] == 0) nxt.push(u);
	}

	assert((int64_t) sorted.size() == n); // true iff input is a dag
	return sorted;
}

void find_path(int64_t src, int64_t dest, std::vector<int64_t> &path, std::vector<std::vector<int64_t>> &dag,
		std::vector<int64_t> &order) {
	std::unordered_map<int64_t, bool> vis;
	std::function<bool(int64_t)> dfs = [&](int64_t current) {
		if (order[current] > order[dest]) return false;
		if (vis[current]) return false;
		path.push_back(current);
		vis[current] = true;
		if (current == dest) return true;
		for (int64_t nxt: dag[current]) {
			if (dfs(nxt)) return true;
		}
		path.pop_back();
		return false;
	};
	assert(dfs(src));
}
