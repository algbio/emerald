#include <string>
#include <iostream>
#include <map>
#include <vector>

#include "draw_subgraph.h"

std::string draw_subgraph(const int64_t IDX, const int64_t n, const int64_t m, const Dag &d, const std::vector<std::vector<mpq_class>> &ratios, mpq_class alpha, const std::string &a, const std::string &b) {
	auto has = [&](int64_t i, int64_t j) {
		// check if pair (i, j) is part of the suboptimal dag
		//return d.trans.find(std::make_pair(i, j)) != d.trans.end();
		if (d.trans.find(std::make_pair(i, j)) != d.trans.end()) {
			std::array<int64_t, 3> def = { -1, -1, -1 };
			std::array<int64_t, 3> our = d.trans.at(std::make_pair(i, j));
			assert(our != def);
			return true;
		}
		return false;
	};
	auto is_safe = [&](int64_t i1, int64_t j1, int64_t i2, int64_t j2) {
		if (!has(i1, j1) || !has(i2, j2)) return false;
		for (int64_t a1: d.trans.at(std::make_pair(i1, j1))) if (a1 > -1) {
			assert(ratios[a1].size() == d.adj[a1].size());
			for (int64_t a2: d.trans.at(std::make_pair(i2, j2))) if (a2 > -1) {
				for (int64_t j = 0; j < (int64_t) d.adj[a1].size(); j++) if (d.adj[a1][j] == a2)
					if (ratios[a1][j] >= alpha) return true;
			}
		}
		return false;
	};
	auto is_opt = [&](int64_t i1, int64_t j1, int64_t i2, int64_t j2) {
		if (!has(i1, j1) || !has(i2, j2)) return false;
		for (int64_t a1: d.trans.at(std::make_pair(i1, j1))) if (a1 > -1) {
			for (int64_t a2: d.trans.at(std::make_pair(i2, j2))) if (a2 > -1) {
				bool is_edge = false;
				for (int64_t nxt: d.adj[a1]) if (nxt == a2) is_edge = true;
				if (!is_edge) continue;
				if (d.in_optimal.at(std::make_pair(a1, a2))) return true;
			}
		}
		return false;
	};
	auto outside = [&](int64_t i, int64_t j) {
		return i < 0 || i >= n || j < 0 || j >= m;
	};
	std::string output = "graph fasta_" + std::to_string(IDX) + " {\nnode [shape=none, width=.1, height=.025, fixedsize=true];\n";//rankdir=LR;\n";
	std::vector<std::string> nodes, edges;
	for (int64_t i = 0; i < n; i++) {
		for (int64_t j = 0; j < m; j++) {
			std::string node_label;
			if (i == 0 && j > 0) node_label = std::string(1, b[j - 1]);
			else if (i > 0 && j == 0) node_label = std::string(1, a[i - 1]);
			else node_label = "";
			//else node_label = std::to_string(i) + "_" + std::to_string(j);
			std::string node = "\"" + std::to_string(i) + "_" + std::to_string(j) + "\" [label=\"" + node_label + "\"";
			if ((i > 0 && j > 0) || (i == 0 && j == 0))
				node += (has(i, j) ? ", style=filled, shape=circle]" : ", style=invis]");
			else
				node += "]";
			nodes.push_back(node);

			std::vector<std::pair<int64_t, int64_t>> from = { {i - 1, j - 1}, {i - 1, j}, {i, j - 1} };
			for (int dir = 0; dir < 3; dir++) {
				auto [l, k] = from[dir];
				if (outside(l, k)) continue;
				bool is_edge = has(i, j) && has(l, k) && d.trans.at(std::make_pair(i, j))[dir] != -1;
				std::string edge = "\"" + std::to_string(l) + "_" + std::to_string(k) + "\" -- \"" + std::to_string(i) + "_" + std::to_string(j) + "\"";
				if (!has(i, j) || !has(l, k) || !is_edge) {
					edge += " [style=invis]";
				} else {
					//int s;
					//if ((s = is_safe(l, k, i, j))) edge += ((s % 2) ? " [color=green, penwidth=5]" : " [color=green, penwidth=5]");
					//else
					if (is_opt(l, k, i, j)) edge += " [penwidth=3]";
					else edge += " [color=grey, penwidth=2]";
				}
				edges.push_back(edge);
			}
		}
	}

	for (int i = 0, k = 0; i < n; i++) {
		output += "subgraph row_" + std::to_string(i) + " {\n";
		output += "rank=same;\n";
		for (int j = 0; j < m; j++) {
			output += nodes[k++] + ";\n";
		}
		output += "}\n";
	}
	for (std::string edge: edges) output += edge + ";\n";
	output += "}\n";

	return output;
}

