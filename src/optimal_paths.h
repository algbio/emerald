#pragma once
#include <vector>
#include <map>
#include <string>
#include <array>
#include <assert.h>

#include <gmpxx.h>


// Dag here is supposed to be the subgraph constructed with the gen_dag function
struct Dag {
	std::vector<std::vector<int64_t>> adj;
	std::vector<std::vector<std::pair<int64_t, int64_t>>> adj_costs; // pair (vertex, cost)
	int64_t src, sink;
	std::map<std::pair<int64_t, int64_t>, std::array<int64_t, 3>> trans;
	std::map<int64_t, std::pair<int64_t, int64_t>> transr;
	std::map<std::pair<int64_t, int64_t>, bool> in_optimal;
};

struct Node {
	int64_t N_index, M_index, type;
	int64_t cost;
	Node (int64_t N_index, int64_t M_index, int64_t type, int64_t cost) : N_index(N_index),
			M_index(M_index), type(type), cost(cost)
	{}
};

// Check whether edge is internal or external
bool is_internal(Dag &d, const int64_t v, const int64_t u);

// Check whether there is a gap for sequence S.
// If for_representative is true, then S is equal to the representative (std::string &a)
// If for_representative is false, then S is equal to the cluster member (std::string &b)
bool is_gap(Dag &d, const int64_t v, const int64_t u, bool for_representative);

// Print suboptimal alignments into a fasta file
std::vector<std::vector<int64_t>> alignments_into_fasta(int64_t print_alignments, Dag &d, const std::string &a, const std::string &b, const std::string &fasta_file);

// Construct alignment paths
std::vector<std::vector<std::vector<std::vector<Node>>>> build_dp_matrix(const std::string &a,
		const std::string &b, const int64_t GAP_COST, const int64_t START_GAP, const int64_t cost_matrix[21][21], int64_t sign);

// Construct optimal alignment score matrix
std::vector<std::vector<std::vector<int64_t>>>
opt_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj, bool dir);

// Returns the score of a randomly chosen alignment
int64_t score_of_random_alignment(const std::vector<std::vector<std::vector<std::vector<Node>>>> &adj);

// Find the sub-graph of the alignment paths with (sub-)optimal paths
template<class K>
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
	auto check_opt = [&](const int64_t k, const int64_t OPT, const int64_t delta) {
		return k == OPT;
	};

	for (int64_t i = 0; i <= n; i++) for (int64_t j = 0; j <= m; j++) for (int64_t k = 2; k >= 0; k--) {
		if (!check_th(dp[i][j][k] + dpr[i][j][k], OPT, delta)) continue;
		add_node(Node(i, j, k, 0));
		for (const Node &nxt: e[i][j][k]) {
			if (!check_th(dp[nxt.N_index][nxt.M_index][nxt.type] + dpr[nxt.N_index][nxt.M_index][nxt.type], OPT, delta)) continue;
			if (check_th(dpr[nxt.N_index][nxt.M_index][nxt.type] + dp[i][j][k] + nxt.cost, OPT, delta)) {
				arcs++;
				add_node(nxt);
				in_optimal[std::make_pair(trans[std::make_pair(i, j)][k],
						trans[std::make_pair(nxt.N_index, nxt.M_index)][nxt.type])] = check_opt(dpr[nxt.N_index][nxt.M_index][nxt.type] + dp[i][j][k] + nxt.cost, OPT, delta);
				adj[trans[std::make_pair(i, j)][k]].push_back(trans[std::make_pair(nxt.N_index, nxt.M_index)][nxt.type]);
				adj_costs[trans[std::make_pair(i, j)][k]].emplace_back((trans[std::make_pair(nxt.N_index, nxt.M_index)][nxt.type]), nxt.cost);
			}
		}
	}

	if (verbose_flag) {
		std::cout << "TOTAL: " << 3*(n+1)*(m+1) << std::endl;
		std::cout << "CURRENT: " << current << std::endl;
		std::cout << "Number of arcs: " << arcs << std::endl;
		std::cout << OPT - delta << std::endl;
		if (subpaths) std::cout << "Found subpaths: " << subpaths << std:: endl;
	}

	return { adj, adj_costs, 0, trans[std::make_pair(n, m)][0], trans, transr, in_optimal };
}
