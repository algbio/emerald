#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <string>
#include <assert.h>
#include <array>
#include <stack>
#include <ctime>

#include <fstream>
#include <iostream>

#include "optimal_paths.h"

// translate fasta file letters to amino acid symbols (see http://www.math.utep.edu/Faculty/mleung/bioinformatics/aacodon.html)
//                        A   B  C  D  E   F  G  H  I   J   K   L   M  N   O   P  Q  R   S   T   U   V   W   X   Y   Z
const int64_t LTA[26] = { 0, 20, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 20 };

bool is_internal(Dag &d, const int64_t v, const int64_t u) {
	return d.transr[v] == d.transr[u];
}

bool is_gap(Dag &d, const int64_t v, const int64_t u, bool for_representative) {
	return !is_internal(d, v, u) && (for_representative ? d.transr[v].first == d.transr[u].first : d.transr[v].second == d.transr[u].second);
}

std::vector<std::vector<int64_t>> alignments_into_fasta(int64_t number_of_paths, Dag &d, const std::string &a, const std::string &b, const std::string &fasta_file) {
	std::string alignment = "";
	std::ofstream fasta_stream;
	fasta_stream.open(fasta_file, std::ofstream::app);
	std::cout << "Print alignments: created file ./" << fasta_file << std::endl;
	int64_t al_idx = 0;
	auto print_alignment = [&](std::string s, int64_t c) {
		fasta_stream << s << '\n';
	};

	// calculate set of paths P as the top suboptimal alignments
	int64_t n = (int64_t) d.adj.size();
	std::vector<std::vector<int64_t>> P;
	std::vector<int64_t> costs_of_final_paths;
	std::vector<int64_t> count(n, 0LL);

	std::vector<std::pair<int64_t, std::vector<int64_t>>> B; // pair (cost, path)
	B.push_back(std::make_pair(0LL, std::vector<int64_t>{d.src}));
	while (!B.empty() && count[d.sink] < number_of_paths) {
		std::pop_heap(B.begin(), B.end());
		auto [cost, Pp] = B.back();
		B.pop_back();
		count[Pp.back()]++;
		if (Pp.back() == d.sink)
			P.push_back(Pp), costs_of_final_paths.push_back(cost);
		if (count[Pp.back()] <= number_of_paths) {
			for (auto [v, nxtc]: d.adj_costs[Pp.back()]) {
				std::vector<int64_t> Pv = Pp;
				Pv.push_back(v);
				B.push_back(std::make_pair(cost + nxtc, Pv));
				std::push_heap(B.begin(), B.end());
			}
		}
	}

	int64_t ci = 0;
	for (const std::vector<int64_t> &path: P) {
		std::string alignment_a = "", alignment_b = "";
		for (int64_t i = 1; i < (int64_t) path.size(); i++) {
			std::pair<int64_t, int64_t> p = d.transr[path[i - 1]];
			std::pair<int64_t, int64_t> n = d.transr[path[i]];
			if (n == p) continue; // node internal edge
			alignment_a += (n.first == p.first + 1 ? a[p.first] : '_');
			alignment_b += (n.second == p.second + 1 ? b[p.second] : '_');
		}
		print_alignment(alignment_a, costs_of_final_paths[ci]);
		print_alignment(alignment_b, costs_of_final_paths[ci++]);
		print_alignment("", 0);
	}

	fasta_stream.close();

	return P; // for unit testing purposes
}

std::vector<std::vector<std::vector<std::vector<Node>>>> build_dp_matrix(const std::string &a,
		const std::string &b, const int64_t GAP_COST, const int64_t START_GAP, const int64_t cost_matrix[21][21],
		int64_t sign) {
	int64_t n = (int64_t) a.size();
	int64_t m = (int64_t) b.size();

	// create affine linear gap cost graph (see README for illustration)
	std::vector<std::vector<std::vector<std::vector<Node>>>> adj(n + 1,
			std::vector<std::vector<std::vector<Node>>>(m + 1,
			std::vector<std::vector<Node>>(3)));
	for (int64_t i = 0; i <= n; i++) for (int64_t j = 0; j <= m; j++) {
		if (i < n) {
			if (a[i] < 'A' || a[i] > 'Z') {
				std::cerr << "ERROR: INVALID CHARACTER in the reference string: " << a[i] << std::endl;
			}
			if (LTA[a[i] - 'A'] == -1) {
				std::cerr << "ERROR: WRONG CHARACTER in the reference string: " << a[i] << std::endl;
			}
		}
		if (j < m) {
			if (b[j] < 'A' || b[j] > 'Z') {
				std::cerr << "ERROR: INVALID CHARACTER in the comparing string: " << b[j] << std::endl;
			}
			if (LTA[b[j] - 'A'] == -1) {
				std::cerr << "ERROR: WRONG CHARACTER in the comparing string: " << b[j] << std::endl;
			}
		}

		if (j + 1 <= m) {
			adj[i][j][0].push_back(Node(i, j + 1, 2, sign * (START_GAP + GAP_COST)));
			adj[i][j][2].push_back(Node(i, j + 1, 2, sign * GAP_COST));
		}

		if (i + 1 <= n) {
			adj[i][j][0].push_back(Node(i + 1, j, 1, sign * (START_GAP + GAP_COST)));
			adj[i][j][1].push_back(Node(i + 1, j, 1, sign * GAP_COST));
		}

		if (i + 1 <= n && j + 1 <= m)
			adj[i][j][0].push_back(Node(i + 1, j + 1, 0,
					sign * cost_matrix[LTA[a[i] - 'A']][LTA[b[j] - 'A']]));

		adj[i][j][1].push_back(Node(i, j, 0, 0));
		adj[i][j][2].push_back(Node(i, j, 0, 0));
	}

	return adj;
}

std::vector<std::vector<std::vector<int64_t>>>
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, bool dir) {
	int64_t n = (int64_t) adj.size() - 1;
	assert(n > 0);
	int64_t m = (int64_t) adj[0].size() - 1;
	assert(m > 0);

	std::vector<std::vector<std::vector<int64_t>>> d(n + 1, std::vector<std::vector<int64_t>>(m + 1, std::vector<int64_t>(3, -(1 << 30))));
	(dir ? d[n][m][0] : d[0][0][0]) = 0;
	auto update_dist = [&](int64_t i, int64_t j, int64_t k) {
		for (const Node &node: adj[i][j][k]) {
			d[node.N_index][node.M_index][node.type] = std::max(d[node.N_index][node.M_index][node.type], d[i][j][k] + node.cost);
		}
	};
	if (!dir) {
		for (int64_t i = 0; i <= n; i++) for (int64_t j = 0; j <= m; j++) for (int64_t k = 2; k >= 0; k--) {
			if (d[i][j][k] <= -(1 << 30)) continue;
			update_dist(i, j, k);
		}
	} else {
		for (int64_t i = n; i >= 0; i--) for (int64_t j = m; j >= 0; j--) for (int64_t k = 0; k <= 2; k++) {
			if (d[i][j][k] <= -(1 << 30)) continue;
			update_dist(i, j, k);
		}
	}
	return d;
}

int64_t score_of_random_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj) {
	int64_t n = (int64_t) adj.size() - 1;
	assert(n > 0);
	int64_t m = (int64_t) adj[0].size() - 1;
	assert(m > 0);

	int64_t score = 0;
	std::tuple<int64_t, int64_t, int64_t> v = std::make_tuple(0, 0, 0);
	std::srand(std::time(nullptr));
	while (v != std::make_tuple(n, m, 0)) {
		auto [i, j, t] = v;
		int64_t idx = rand() % adj[i][j][t].size();
		Node nxt = adj[i][j][t][idx];
		score += nxt.cost;
		v = std::make_tuple(nxt.N_index, nxt.M_index, nxt.type);
	}
	return score;
}

Dag gen_dag(const std::string &a, const std::string &b, const int64_t cost_matrix[21][21],
		const int64_t delta, const int64_t GAP_COST, const int64_t START_GAP, bool random_alignment_as_optimal, int &verbose_flag) {
	std::vector<std::vector<std::vector<std::vector<Node>>>> e = build_dp_matrix(a, b, GAP_COST, START_GAP, cost_matrix, 1);
	int64_t n = (int64_t) e.size() - 1;
	assert(n > 0);
	int64_t m = (int64_t) e[0].size() - 1;
	assert(m > 0);

	// find highest scores for each node
	std::vector<std::vector<std::vector<std::vector<Node>>>> er(n + 1, std::vector<std::vector<std::vector<Node>>>(m + 1, std::vector<std::vector<Node>>(3)));
	for (int64_t i = 0; i <= n; i++) for (int64_t j = 0; j <= m; j++) for (int64_t k = 2; k >= 0; k--) {
		for (Node nxt: e[i][j][k])
			er[nxt.N_index][nxt.M_index][nxt.type].push_back(Node(i, j, k, nxt.cost));
	}

	std::vector<std::vector<std::vector<int64_t>>> dp = opt_alignment(e, 0);
	std::vector<std::vector<std::vector<int64_t>>> dpr = opt_alignment(er, 1);
	assert(dp[n][m][0] == dpr[0][0][0]);

	const int64_t OPT = (random_alignment_as_optimal ? score_of_random_alignment(e) : dp[n][m][0]);

	int64_t current = 0;
	std::vector<std::vector<int64_t>> adj(1);
	std::vector<std::vector<std::pair<int64_t, int64_t>>> adj_costs(1);
	std::map<std::pair<int64_t, int64_t>, std::array<int64_t, 3>> trans; // translate to index
	std::map<int64_t, std::pair<int64_t, int64_t>> transr; // translate index to pair
	trans[std::make_pair(0, 0)] = {0, -1, -1};
	transr[0] = std::make_pair(0, 0);
	if (verbose_flag)
		std::cout << "OPT: " << OPT << std::endl;

	auto add_node = [&](const Node &node) {
		if (trans.find(std::make_pair(node.N_index, node.M_index)) == trans.end()) {
			trans[std::make_pair(node.N_index, node.M_index)] = {-1, -1, -1};
		}
		if (trans[std::make_pair(node.N_index, node.M_index)][node.type] == -1) {
			trans[std::make_pair(node.N_index, node.M_index)][node.type] = ++current;
			transr[current] = std::make_pair(node.N_index, node.M_index);
			adj.push_back(std::vector<int64_t>());
			adj_costs.push_back(std::vector<std::pair<int64_t, int64_t>>());
		}
	};


	int64_t subpaths = 0, arcs = 0;
	std::map<std::pair<int64_t, int64_t>, bool> in_optimal;
	auto check_th = [&](const int64_t k, const int64_t OPT, const int64_t delta) {
		if (k < OPT && k >= OPT - delta) subpaths++;
		return k >= OPT - delta && k <= OPT + delta;
	};
}
