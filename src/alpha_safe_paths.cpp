#include <iostream>
#include <functional>
#include <vector>
#include <set>
#include <stack>
#include <unordered_map>
#include <assert.h>

#include <gmpxx.h>

#include "alpha_safe_paths.h"
#include "optimal_paths.h"


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

std::vector<mpz_class> number_of_paths(std::vector<std::vector<int64_t>> &dag) {
	int64_t n = (int64_t) dag.size();
	std::vector<int64_t> sorted = topsort(dag);
	int64_t sink = sorted.back();

	std::vector<mpz_class> am(n, mpz_class(0));
	am[sink] = 1;
	for (int64_t i = n - 1; i >= 0; i--) {
		for (int64_t v: dag[sorted[i]]) am[sorted[i]] += am[v];
	}
	return am;
}

std::vector<std::vector<mpq_class>> path_ratios(Dag &d, std::vector<mpz_class> &am, std::vector<mpz_class> &ram) {
	std::vector<std::vector<int64_t>> &dag = d.adj;
	int64_t n = (int64_t) dag.size();

	for (int64_t i = 0; i < n; i++) assert(am[i] <= am[d.src]);

	std::vector<std::vector<mpq_class>> ratios(n);

	int64_t SRC = d.src;

	for (int64_t i = 0; i < n; i++) for (int64_t v: dag[i]) {
		mpq_class nxt = am[v] * ram[i];
		nxt /= am[SRC];
		ratios[i].push_back(nxt);
	}
	return ratios;
}

std::vector<int64_t> find_alpha_path(Dag &d,
		std::vector<std::vector<mpq_class>> &ratios, mpq_class alpha, bool verbose_flag) {
	std::vector<std::vector<int64_t>> &dag = d.adj;
	int64_t n = (int64_t) ratios.size();
	std::vector<int64_t> sorted = topsort(dag);

	std::vector<std::pair<int64_t, int64_t>> needed;
	for (int64_t i = 0; i < n; i++) {
		int64_t u = sorted[i];
		int64_t am = 0;
		for (int64_t j = 0; j < (int64_t) dag[u].size(); j++) {
			int64_t v = dag[u][j];
			mpq_class d = ratios[u][j];
			assert(d <= 1);
			if (d >= alpha) needed.emplace_back(u, v), am++;
		}
		assert(am <= 1);
	}
	if (verbose_flag) {
		int64_t sed = 0;
		for (auto [u, v]: needed) sed += !is_internal(d, u, v);
		std::cout << "Number of edges in SW: " << (int64_t) needed.size() << '\n';
		std::cout << "Number of safe external edges: " << sed << '\n';
	}

	std::vector<int64_t> order(n);
	for (int64_t i = 0; i < n; i++) order[sorted[i]] = i;

	std::vector<int64_t> path;
	int64_t last = d.src;
	for (auto [u,v]: needed) {
		find_path(last, u, path, dag, order);
		last = v;
	}
	find_path(last, d.sink, path, dag, order);

	return path;
}

