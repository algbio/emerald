#include "wasm_interface.h"
#include "cost_matrices.h"
#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"
#include "write_json.h"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <sstream>

// Helper function to get a predefined cost matrix
const int64_t (*get_cost_matrix(CostMatrixType matrix_type))[21] {
    switch (matrix_type) {
        case BLOSUM62:
            return BLOSUM62_MATRIX;
        case PAM250:
            return PAM250_MATRIX;
        case IDENTITY:
            return IDENTITY_MATRIX;
        default:
            throw std::invalid_argument("Unknown cost matrix type");
    }
}

// Copy matrix values from a vector representation to the 2D array format used by the alignment functions
void copy_matrix_from_vector(const std::vector<std::vector<int64_t>>& src, int64_t dest[21][21]) {
    // First check that the input matrix has the correct size
    if (src.size() != 21) {
        throw std::invalid_argument("Custom matrix must have 21 rows");
    }
    
    for (size_t i = 0; i < src.size(); i++) {
        if (src[i].size() != 21) {
            throw std::invalid_argument("Each row in custom matrix must have 21 columns");
        }
        
        for (size_t j = 0; j < src[i].size(); j++) {
            dest[i][j] = src[i][j];
        }
    }
}

// Implementation for the function with predefined matrix selection
std::string generate_alignment_json_with_matrix(const std::string& representative_sequence,
                                              const std::string& representative_descriptor,
                                              const std::string& member_sequence,
                                              const std::string& member_descriptor,
                                              CostMatrixType matrix_type,
                                              float alpha,
                                              int64_t delta,
                                              int64_t gap_cost,
                                              int64_t start_gap) {
    
    // Create protein objects
    Protein ref(std::string(">" + representative_descriptor));
    ref.sequence = representative_sequence;
    
    Protein mem(std::string(">" + member_descriptor));
    mem.sequence = member_sequence;
    
    // Get the selected cost matrix
    const int64_t (*selected_matrix)[21] = get_cost_matrix(matrix_type);
    
    // Copy to a mutable array for use in the alignment function
    int64_t cost_matrix[21][21];
    for (int i = 0; i < 21; i++) {
        for (int j = 0; j < 21; j++) {
            cost_matrix[i][j] = selected_matrix[i][j];
        }
    }
    
    // Use defaults if not specified
    if (gap_cost == 0) gap_cost = -1;  
    if (start_gap == 0) start_gap = -11;
    
    // Create suboptimal space
    bool random_alignment_as_optimal = false;
    int verbose_flag = 0;
    Dag d = gen_dag(ref.sequence, mem.sequence, cost_matrix, delta, gap_cost, start_gap, random_alignment_as_optimal, verbose_flag);
    
    // Find the number of v-t paths (am) and the number of s-v paths (ram) for all nodes v
    std::vector<std::vector<int64_t>> adj = d.adj;
    std::vector<std::vector<int64_t>> radj((int64_t) adj.size());
    for (int64_t i = 0; i < (int64_t) adj.size(); i++) {
        for (int64_t v: adj[i]) {
            radj[v].push_back(i);
        }
    }
    
    std::vector<mpz_class> am = number_of_paths(adj);
    std::vector<mpz_class> ram = number_of_paths(radj);

    // Find path ratios
    std::vector<std::vector<mpq_class>> ratios = path_ratios(d, am, ram);
    
    // Find a path that contains all safety windows
    std::vector<int64_t> path = find_alpha_path(d, ratios, alpha, verbose_flag);
    
    // Find safety windows
    std::vector<mpq_class> r = find_ratios(path, adj, ratios);
    auto [swindows, window_ratios, number_of_edges] = safety_windows(am, ram, path, alpha);
    
    // Extract coordinate windows
    std::vector<std::pair<int64_t, int64_t>> windows, windowsp;
    for (int64_t i = 0; i < (int64_t) swindows.size(); i++) {
        auto [LT, RT] = swindows[i];
        int64_t L = d.transr[LT].first, R = d.transr[RT].first;
        int64_t Lp = d.transr[LT].second, Rp = d.transr[RT].second;

        windows.emplace_back(L, R);
        windowsp.emplace_back(Lp, Rp);
    }
    
    // Generate alignment strings
    auto [alignment_ref, alignment_mem] = generate_alignment_strings(d, path, ref, mem);
    
    // Generate JSON
    std::stringstream json_stream;
    write_json_to_stream(json_stream, ref, mem, d, windows, windowsp, ratios, alignment_ref, alignment_mem);
    return json_stream.str();
}

// Implementation for the function with custom cost matrix
std::string generate_alignment_json_with_custom_matrix(const std::string& representative_sequence,
                                                     const std::string& representative_descriptor,
                                                     const std::string& member_sequence,
                                                     const std::string& member_descriptor,
                                                     const std::vector<std::vector<int64_t>>& custom_matrix,
                                                     float alpha,
                                                     int64_t delta,
                                                     int64_t gap_cost,
                                                     int64_t start_gap) {
    
    // Create protein objects
    Protein ref(std::string(">" + representative_descriptor));
    ref.sequence = representative_sequence;
    
    Protein mem(std::string(">" + member_descriptor));
    mem.sequence = member_sequence;
    
    // Copy custom matrix to the required 21x21 format
    int64_t cost_matrix[21][21];
    copy_matrix_from_vector(custom_matrix, cost_matrix);
    
    // Use defaults if not specified
    if (gap_cost == 0) gap_cost = -1;  
    if (start_gap == 0) start_gap = -11;
    
    // Create suboptimal space
    bool random_alignment_as_optimal = false;
    int verbose_flag = 0;
    Dag d = gen_dag(ref.sequence, mem.sequence, cost_matrix, delta, gap_cost, start_gap, random_alignment_as_optimal, verbose_flag);
    
    // Find the number of v-t paths (am) and the number of s-v paths (ram) for all nodes v
    std::vector<std::vector<int64_t>> adj = d.adj;
    std::vector<std::vector<int64_t>> radj((int64_t) adj.size());
    for (int64_t i = 0; i < (int64_t) adj.size(); i++) {
        for (int64_t v: adj[i]) {
            radj[v].push_back(i);
        }
    }
    
    std::vector<mpz_class> am = number_of_paths(adj);
    std::vector<mpz_class> ram = number_of_paths(radj);

    // Find path ratios
    std::vector<std::vector<mpq_class>> ratios = path_ratios(d, am, ram);
    
    // Find a path that contains all safety windows
    std::vector<int64_t> path = find_alpha_path(d, ratios, alpha, verbose_flag);
    
    // Find safety windows
    std::vector<mpq_class> r = find_ratios(path, adj, ratios);
    auto [swindows, window_ratios, number_of_edges] = safety_windows(am, ram, path, alpha);
    
    // Extract coordinate windows
    std::vector<std::pair<int64_t, int64_t>> windows, windowsp;
    for (int64_t i = 0; i < (int64_t) swindows.size(); i++) {
        auto [LT, RT] = swindows[i];
        int64_t L = d.transr[LT].first, R = d.transr[RT].first;
        int64_t Lp = d.transr[LT].second, Rp = d.transr[RT].second;

        windows.emplace_back(L, R);
        windowsp.emplace_back(Lp, Rp);
    }
    
    // Generate alignment strings
    auto [alignment_ref, alignment_mem] = generate_alignment_strings(d, path, ref, mem);
    
    // Generate JSON
    std::stringstream json_stream;
    write_json_to_stream(json_stream, ref, mem, d, windows, windowsp, ratios, alignment_ref, alignment_mem);
    return json_stream.str();
}
