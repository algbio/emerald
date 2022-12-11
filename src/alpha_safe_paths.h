#pragma once
#include <vector>
#include <gmpxx.h>

#include "optimal_paths.h"

// given a dag of optimal paths, find a path with almost safe (>= alpha) paths

std::vector<int64_t> topsort(std::vector<std::vector<int64_t>> &dag);

// Find path from src to dest. Assertion error in case if such path does not exist
void find_path(int64_t src, int64_t dest, std::vector<int64_t> &path,
		std::vector<std::vector<int64_t>> &dag, std::vector<int64_t> &order);

// For each vertex, save the amount of paths starting from the vertex to the sink
template<class T>
std::vector<T> number_of_paths(std::vector<std::vector<int64_t>> &dag) {
	int64_t n = (int64_t) dag.size();
	std::vector<int64_t> sorted = topsort(dag);
	int64_t sink = sorted.back();

	//std::vector<mpz_class> am(n, mpz_class(0));
	std::vector<T> am(n, T(0));
	am[sink] = 1;
	for (int64_t i = n - 1; i >= 0; i--) {
		for (int64_t v: dag[sorted[i]]) am[sorted[i]] += am[v];
	}
	return am;
}

// For each edge, calculate the % of s-t paths they are part in
template<class T, class K>
std::vector<std::vector<K>> path_ratios(Dag &d, std::vector<T> &am, std::vector<T> &ram) {
	std::vector<std::vector<int64_t>> &dag = d.adj;
	int64_t n = (int64_t) dag.size();

	for (int64_t i = 0; i < n; i++) assert(am[i] <= am[d.src]);

	std::vector<std::vector<K>> ratios(n);

	int64_t SRC = d.src;

	for (int64_t i = 0; i < n; i++) for (int64_t v: dag[i]) {
		// nxt = am[v] * ram[i] / am[SRC]
		//mpq_class nxt = am[v] * ram[i];
		K nxt = am[v] * ram[i];
		nxt /= am[SRC];
		ratios[i].push_back(nxt);
	}
	return ratios;
}

// Find s--t path that contains all edges with occurence ratio >= alpha. Might fail if alpha <= 0.5
template<class K>
std::vector<int64_t> find_alpha_path(Dag &d,
		std::vector<std::vector<K>> &ratios, K alpha, bool verbose_flag) {
	std::vector<std::vector<int64_t>> &dag = d.adj;
	int64_t n = (int64_t) ratios.size();
	std::vector<int64_t> sorted = topsort(dag);

	std::vector<std::pair<int64_t, int64_t>> needed;
	for (int64_t i = 0; i < n; i++) {
		int64_t u = sorted[i];
		int64_t am = 0;
		for (int64_t j = 0; j < (int64_t) dag[u].size(); j++) {
			int64_t v = dag[u][j];
			K d = ratios[u][j];
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

