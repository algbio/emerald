#include "safety_window_calculator.h"
#include <iostream>
#include <queue>
#include <functional>
#include <algorithm>
#include <assert.h>

// Helper for using pairs as map keys
struct PairHash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

std::vector<SafetyWindow> SafetyWindowCalculator::computeSafetyWindows(
    const AlignmentResult& alignment_result,
    const DPMatrix::ScoreMatrix& forward_scores,
    const DPMatrix::ScoreMatrix& backward_scores,
    size_t n, size_t m) {
    
    if (alignment_result.delta_neighborhood.empty()) {
        if (config_.verbose) {
            std::cout << "No edges in delta neighborhood, no safety windows." << std::endl;
        }
        return {};
    }
    if (config_.verbose) {
        std::cout << "Computing safety windows from delta neighborhood with "
              << alignment_result.delta_neighborhood.size() << " edges." << std::endl;
    }

    // Build a graph representation from the delta neighborhood
    std::unordered_map<size_t, std::vector<size_t>> graph, reverse_graph;
    std::vector<size_t> topo_order;
    std::vector<size_t> alignment_with_all_safe_edges;
    std::unordered_map<std::pair<size_t, size_t>, int64_t, PairHash> edge_cost;
    std::unordered_map<size_t, std::pair<size_t, size_t>> node_position;
    
    // Function to convert (i,j,state) to unique node ID
    auto nodeToId = [m](size_t i, size_t j, int state) -> size_t {
        return i * (m+1) * 3 + j * 3 + state;
    };
    
    // Function to extract (i,j,state) from node ID
    auto idToNode = [m](size_t id) -> std::tuple<size_t, size_t, int> {
        int state = id % 3;
        size_t j = ((id - state) / 3) % (m+1);
        size_t i = (id - state - 3*j) / ((m+1) * 3);
        return {i, j, state};
    };
    
    // Process all edges in delta neighborhood to build graph
    for (const auto& edge : alignment_result.delta_neighborhood) {
        size_t src_id = nodeToId(
            edge.source.getI(), 
            edge.source.getJ(), 
            static_cast<int>(edge.source.getState())
        );
        
        size_t dst_id = nodeToId(
            edge.target.getI(), 
            edge.target.getJ(), 
            static_cast<int>(edge.target.getState())
        );
        
        // Store graph structure
        graph[src_id].push_back(dst_id);
        reverse_graph[dst_id].push_back(src_id);
        
        // Store edge cost and node positions
        edge_cost[{src_id, dst_id}] = edge.cost;
        node_position[src_id] = {edge.source.getI(), edge.source.getJ()};
        node_position[dst_id] = {edge.target.getI(), edge.target.getJ()};
    }

    // Calculate paths through each edge using dynamic programming
    size_t source = nodeToId(0, 0, static_cast<int>(AlignmentState::MATCH));
    size_t sink = nodeToId(n, m, static_cast<int>(AlignmentState::MATCH));
    
    // Step 1: Count paths from source to each node
    std::unordered_map<size_t, double> paths_to;
    paths_to[source] = 1.0;
    
    // Topological sort
    std::unordered_set<size_t> visited;
    std::function<void(size_t)> dfs = [&](size_t node) {
        visited.insert(node);
        if (graph.count(node)) {
            for (size_t next : graph[node]) {
                if (!visited.count(next)) {
                    dfs(next);
                }
            }
        }
        topo_order.push_back(node);
    };
    
    dfs(source);
    std::reverse(topo_order.begin(), topo_order.end());

    // Build a rank map from the forward topological order (before clearing it later)
    std::unordered_map<size_t, size_t> topo_rank;
    for (size_t idx = 0; idx < topo_order.size(); ++idx) {
        topo_rank[topo_order[idx]] = idx;
    }

    // Compute paths to each node using topological ordering
    for (size_t node : topo_order) {
        assert (reverse_graph.count(node) || node == source);
        if (reverse_graph.count(node)) {
            for (size_t prev : reverse_graph[node]) {
                paths_to[node] += paths_to[prev];
            }
        }
    }
    
    // Step 2: Count paths from each node to sink
    std::unordered_map<size_t, double> paths_from;
    paths_from[sink] = 1.0;
    
    for (auto it = topo_order.rbegin(); it != topo_order.rend(); ++it) {
        size_t node = *it;
        assert (graph.count(node) || node == sink);
        if (graph.count(node)) {
            for (size_t next : graph[node]) {
                paths_from[node] += paths_from[next];
            }
        }
    }
    
    // Total number of paths from source to sink
    double total_paths = paths_to[sink];
    
    if (total_paths < 1e-10) {
        // No valid paths or numerical issues
        return {};
    }
    
    if (config_.verbose) {
        std::cout << "Total paths from source to sink: " << total_paths << std::endl;
    }
    
    // Step 3: Find edges that appear in at least alpha fraction of all paths
    double alpha = config_.alpha;
    std::vector<std::pair<size_t, size_t>> high_freq_edges;

    for (const auto& [src, dests] : graph) {
        for (size_t dst : dests) {
            // Edge frequency = (paths to src) * (paths from dst) / (total paths)
            double freq = (paths_to[src] * paths_from[dst]) / total_paths;

            if (freq >= alpha) {
                high_freq_edges.emplace_back(src, dst);
            }
        }
    }

    if (config_.verbose) {
        std::cout << "Found " << high_freq_edges.size() 
                  << " high-frequency edges with alpha=" << alpha << std::endl;
    }

    // Sort edges by DAG order (source node’s topo rank), then by (i_src, j_src) for stability
    std::sort(high_freq_edges.begin(), high_freq_edges.end(), [&](const auto& a, const auto& b) {
        const auto& [u, v] = a;
        const auto& [x, y] = b;
        if (u != x) return topo_rank[u] < topo_rank[x];
        return topo_rank[v] < topo_rank[y];
    });

    // Step 4: Group high-frequency edges into windows, iterating in DAG order
    std::vector<SafetyWindow> windows;
    if (high_freq_edges.empty()) {
        return windows;
    }
    auto outside = [&](const size_t &L, const size_t &R) {
        if (windows.empty()) return false;
        auto [bL, bR, _] = windows.back();
        return topo_rank[bL] >= topo_rank[L] && topo_rank[bR] <= topo_rank[R];
    };
    auto inside = [&](const size_t &L, const size_t &R) {
        if (windows.empty()) return false;
        auto [bL, bR, _] = windows.back();
        return topo_rank[L] >= topo_rank[bL] && topo_rank[R] <= topo_rank[bR];
    };
    double a = (paths_to[high_freq_edges[0].first] * paths_from[high_freq_edges[0].second]) / total_paths;
	for (size_t L = 0, R = 0; R < high_freq_edges.size();
            a = a * (R + 1 == high_freq_edges.size() ? 1 : paths_from[high_freq_edges[R + 1].second]) / paths_from[high_freq_edges[R].second], R++) {
		while (L < R && a < alpha) {
			a = a * paths_to[high_freq_edges[L + 1].first] / paths_to[high_freq_edges[L].first];
			L++;
		}
		while (outside(high_freq_edges[L].first, high_freq_edges[R].second)) windows.pop_back();
		if (L < R && !inside(high_freq_edges[L].first, high_freq_edges[R].second)) windows.emplace_back(high_freq_edges[L].first, high_freq_edges[R].second, a);
	}  
    
    if (config_.verbose) {
        std::cout << std::endl;
        for (size_t i = 0; i < high_freq_edges.size(); ++i) {
            const auto& [src, dst] = high_freq_edges[i];
            const auto& [src_i, src_j, src_state] = idToNode(src);
            const auto& [dst_i, dst_j, dst_state] = idToNode(dst);
            double freq = (paths_to[src] * paths_from[dst]) / total_paths;
            std::cout << "High-freq edge " << i+1 << ": " << src << '-' << dst << " "
                      << "(" << src_i << "," << src_j << "," << src_state << ") -> "
                      << "(" << dst_i << "," << dst_j << "," << dst_state << ") "
                      << "Freq: " << freq * 100 << '%' << std::endl;
        }
        std::cout << std::endl;
        std::cout << "Found " << windows.size() << " safety windows with alpha=" 
                  << alpha << std::endl;
        for (const auto& window : windows) {
            const auto& [rep_start, mem_start, _] = idToNode(window.start_pos);
            const auto& [rep_end, mem_end, __] = idToNode(window.end_pos);
            std::cout << "Window graph indices: " << window.start_pos << "-" << window.end_pos
                      << ", Representative: " << rep_start << "-" << rep_end
                      << ", Member: " << mem_start << "-" << mem_end
                      << ", Safety: " << window.ratio * 100 << '%' << std::endl;
        }
    }
    
    return windows;
}