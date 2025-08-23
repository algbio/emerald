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
        return {};
    }

    // Build a graph representation from the delta neighborhood
    std::unordered_map<size_t, std::vector<size_t>> graph;
    std::unordered_map<size_t, std::vector<size_t>> reverse_graph;
    std::unordered_map<std::pair<size_t, size_t>, int64_t, PairHash> edge_cost;
    std::unordered_map<size_t, std::pair<size_t, size_t>> node_position;
    
    // Function to convert (i,j,state) to unique node ID
    auto nodeToId = [m](size_t i, size_t j, int state) {
        return i * (m+1) * 3 + j * 3 + state;
    };
    
    // Function to extract (i,j,state) from node ID
    auto idToNode = [m](size_t id) -> std::tuple<size_t, size_t, int> {
        int state = id % 3;
        size_t j = (id / 3) % (m+1);
        size_t i = id / ((m+1) * 3);
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
    std::vector<size_t> topo_order;
    std::unordered_set<size_t> visited;
    
    // DFS for topological sort
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
    
    topo_order.clear();
    visited.clear();
    
    // Reverse DFS for reverse topological sort
    std::function<void(size_t)> reverse_dfs = [&](size_t node) {
        visited.insert(node);
        if (reverse_graph.count(node)) {
            for (size_t prev : reverse_graph[node]) {
                if (!visited.count(prev)) {
                    reverse_dfs(prev);
                }
            }
        }
        topo_order.push_back(node);
    };
    
    reverse_dfs(sink);
    std::reverse(topo_order.begin(), topo_order.end());
    
    // Compute paths from each node using reverse topological ordering
    for (size_t node : topo_order) {
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
    std::vector<std::tuple<size_t, size_t, size_t, size_t, double>> high_freq_edges;
    std::vector<std::pair<size_t, size_t>> order; // (topo_rank[src], index into high_freq_edges)

    for (const auto& [src, dests] : graph) {
        for (size_t dst : dests) {
            // Edge frequency = (paths to src) * (paths from dst) / (total paths)
            double freq = (paths_to[src] * paths_from[dst]) / total_paths;

            if (freq >= alpha) {
                auto [i_src, j_src, _] = idToNode(src);
                auto [i_dst, j_dst, __] = idToNode(dst);
                high_freq_edges.emplace_back(i_src, j_src, i_dst, j_dst, freq);
                order.emplace_back(topo_rank[src], high_freq_edges.size() - 1);
            }
        }
    }

    // Sort edges by DAG order (source node’s topo rank), then by (i_src, j_src) for stability
    std::sort(order.begin(), order.end(), [&](const auto& a, const auto& b) {
        if (a.first != b.first) return a.first < b.first;
        const auto& ea = high_freq_edges[a.second];
        const auto& eb = high_freq_edges[b.second];
        if (std::get<0>(ea) != std::get<0>(eb)) return std::get<0>(ea) < std::get<0>(eb); // i_src
        return std::get<1>(ea) < std::get<1>(eb); // j_src
    });

    // Step 4: Group high-frequency edges into windows, iterating in DAG order
    std::vector<SafetyWindow> windows;
    if (order.empty()) {
        return windows;
    }

    auto [start_i, start_j, end_i, end_j, sum_conf] = high_freq_edges[order[0].second];
    int edge_count = 1;

    size_t min_window_size = config_.min_window_size;
    size_t window_gap = config_.safety_window_size;

    for (size_t k = 1; k < order.size(); ++k) {
        const auto& e = high_freq_edges[order[k].second];
        size_t i_src = std::get<0>(e);
        size_t j_src = std::get<1>(e);
        size_t i_dst = std::get<2>(e);
        size_t j_dst = std::get<3>(e);
        double freq   = std::get<4>(e);

        // If this edge continues the current window (within gap tolerance)
        if (i_src <= end_i + window_gap) {
            end_i = std::max(end_i, i_dst);
            end_j = std::max(end_j, j_dst);
            sum_conf += freq;
            edge_count++;
        } else {
            if (end_i - start_i + 1 >= min_window_size) {
                windows.emplace_back(start_i, end_i, start_j, end_j, sum_conf / edge_count);
            }
            start_i = i_src; start_j = j_src;
            end_i = i_dst;   end_j = j_dst;
            sum_conf = freq; edge_count = 1;
        }
    }

    if (end_i - start_i + 1 >= min_window_size) {
        windows.emplace_back(start_i, end_i, start_j, end_j, sum_conf / edge_count);
    }
    
    if (config_.verbose) {
        std::cout << "Found " << windows.size() << " safety windows with alpha=" 
                  << alpha << std::endl;
        for (const auto& window : windows) {
            std::cout << "Window: " << window.start_pos << "-" << window.end_pos
                      << " (target: " << window.target_start << "-" << window.target_end
                      << "), confidence: " << window.confidence << std::endl;
        }
    }
    
    return windows;
}