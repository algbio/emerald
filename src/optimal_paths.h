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
Dag gen_dag(const std::string &a, const std::string &b, const int64_t cost_matrix[21][21],
		const int64_t delta, const int64_t GAP_COST, const int64_t START_GAP, bool random_alignment_as_optimal, int &verbose_flag);
