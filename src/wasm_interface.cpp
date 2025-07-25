#include <algorithm>
#include <vector>
#include <string>
#include <cstring>
#include <gmpxx.h>
#include <map>
#include <unordered_map>
#include <sstream>
#include <fstream>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>  // This is the crucial include for EMSCRIPTEN_BINDINGS
#else
#define EMSCRIPTEN_KEEPALIVE
#endif

#include "wasm_interface.h"
#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"
#include "write_json.h"

// Helper function to generate alignment strings from a path
std::pair<std::string, std::string> generate_alignment_strings(Dag &d, const std::vector<int64_t> &path, const Protein &a, const Protein &b) {
    std::string alignment_a = "", alignment_b = "";
    for (size_t i = 1; i < path.size(); i++) {
        std::pair<int64_t, int64_t> p = d.transr[path[i - 1]];
        std::pair<int64_t, int64_t> n = d.transr[path[i]];
        if (n == p) continue; // node internal edge
        alignment_a += (n.first == p.first + 1 ? a.sequence[p.first] : '-');
        alignment_b += (n.second == p.second + 1 ? b.sequence[p.second] : '-');
    }
    return std::make_pair(alignment_a, alignment_b);
}

// Implementation of the original C++ style function
std::string generate_alignment_json(const std::string& representative_sequence,
                                  const std::string& representative_descriptor,
                                  const std::string& member_sequence,
                                  const std::string& member_descriptor,
                                  float alpha,
                                  int64_t delta,
                                  int64_t gap_cost,
                                  int64_t start_gap) {
    
    // Create protein objects
    Protein ref(std::string(">" + representative_descriptor));
    ref.sequence = representative_sequence;
    
    Protein mem(std::string(">" + member_descriptor));
    mem.sequence = member_sequence;
    
    // Default cost matrix (BLOSUM62)
    int64_t SP = -1;
    int64_t cost_matrix[21][21] = {
        // Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val  Def
        {   4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0,  SP },
        {  -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3,  SP },
        {  -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3,  SP },
        {  -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3,  SP },
        {   0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1,  SP },
        {  -1,   1,   0,   0,  -3,   5,   2,  -2,   0,  -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2,  SP },
        {  -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2,  SP },
        {   0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3,  SP },
        {  -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3,  SP },
        {  -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,  -0,  -3,  -2,  -1,  -3,  -1,   3,  SP },
        {  -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4,  -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1,  SP },
        {  -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,  SP },
        {  -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1,  SP },
        {  -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1,  SP },
        {  -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2,  SP },
        {   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2,  SP },
        {   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0,  SP },
        {  -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2,  11,   2,  -3,  SP },
        {  -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1,  SP },
        {   0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4,  SP },
        {  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP },
    };
    
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

// C-style wrapper implementation
extern "C" {
const char* generate_alignment_json_c(const char* refSeq, const char* refDesc, 
                                    const char* memSeq, const char* memDesc,
                                    float alpha, int64_t delta, int64_t gapCost, int64_t startGap) {
    static std::string result_str;
    
    // Call the C++ version
    result_str = generate_alignment_json(refSeq, refDesc, memSeq, memDesc, 
                                      alpha, delta, gapCost, startGap);
    
    // Return pointer to the static buffer (will persist until next call)
    return result_str.c_str();
}
}

#ifdef __EMSCRIPTEN__
// Embind declarations
EMSCRIPTEN_BINDINGS(emerald_module) {
    using namespace emscripten;
    
    // Original direct JSON return function
    function("generateAlignmentJson", &generate_alignment_json);
}
#endif

// Standalone function that doesn't rely on main program's argument parsing
std::string generateAlignmentJson(const std::string& refSeq, 
                                const std::string& refDesc,
                                const std::string& memSeq, 
                                const std::string& memDesc,
                                float alpha = 0.75, 
                                int64_t delta = 0,
                                int64_t gapCost = -1, 
                                int64_t startGap = -11) {
    // This function is just a wrapper around the main generate_alignment_json function
    return generate_alignment_json(refSeq, refDesc, memSeq, memDesc, alpha, delta, gapCost, startGap);
}
