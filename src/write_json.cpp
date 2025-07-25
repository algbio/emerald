#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <gmpxx.h>
#include <map>  // Add this explicit include

// Include headers that define Protein and Dag structures
#include "alpha_safe_paths.h"  // Assuming this contains the Dag definition
#include "optimal_paths.h"     // Assuming this or another header defines Protein

#include "write_json.h"

void write_json_to_stream(std::ostream &output_stream, Protein &ref, Protein &mem, Dag &d, 
                         std::vector<std::pair<int64_t, int64_t>> &windows, 
                         std::vector<std::pair<int64_t, int64_t>> &windowsp, 
                         std::vector<std::vector<mpq_class>> &ratio,
                         const std::string &alignment_ref,
                         const std::string &alignment_mem)
{
    output_stream << "{\n";
    {
        output_stream << "\t\"representative_descriptor\": \"" << ref.descriptor << "\",\n";
        output_stream << "\t\"representative_string\": \"" << ref.sequence << "\",\n";
        output_stream << "\t\"mem_descriptor\": \"" << mem.descriptor << "\",\n";
        output_stream << "\t\"member_string\": \"" << mem.sequence << "\",\n";
        output_stream << "\t\"alignment_representative\": \"" << alignment_ref << "\",\n";
        output_stream << "\t\"alignment_member\": \"" << alignment_mem << "\",\n";

        // Group edges by coordinates instead of node IDs
        std::map<std::pair<int64_t, int64_t>, std::map<std::pair<int64_t, int64_t>, double>> coord_graph;
        
        for (int64_t i = 0; i < (int64_t)d.adj.size(); i++) {
            auto coords = d.transr.at(i);
            
            for (int64_t j = 0; j < (int64_t)d.adj[i].size(); j++) {
                auto adj_coords = d.transr.at(d.adj[i][j]);
                double edge_ratio = ratio[i][j].get_d();
                
                // Sum ratios if multiple internal edges connect the same coordinate pairs
                coord_graph[coords][adj_coords] += edge_ratio;
            }
        }

        output_stream << "\t\"alignment_graph\": [\n";
        {
            bool first = true;
            for (const auto &[from_coords, edges] : coord_graph) {
                if (!first) output_stream << ",\n";
                first = false;
                
                output_stream << "\t\t{\n";
                output_stream << "\t\t\t\"from\": [" << from_coords.first << ", " << from_coords.second << "],\n";
                output_stream << "\t\t\t\"edges\": [";
                
                bool first_edge = true;
                for (const auto &[to_coords, total_ratio] : edges) {
                    if (!first_edge) output_stream << ", ";
                    first_edge = false;
                    output_stream << "[" << to_coords.first << ", " << to_coords.second << ", " << total_ratio << "]";
                }
                output_stream << "]\n\t\t}";
            }
            output_stream << "\n\t],\n";
        }

        output_stream << "\t\"windows_representative\": [";
        bool first_window = true;
        for (auto [l, r] : windows) {
            if (!first_window) output_stream << ", ";
            first_window = false;
            output_stream << "[" << l << ", " << r << "]";
        }
        output_stream << "],\n";
        
        output_stream << "\t\"windows_member\": [";
        first_window = true;
        for (auto [l, r] : windowsp) {
            if (!first_window) output_stream << ", ";
            first_window = false;
            output_stream << "[" << l << ", " << r << "]";
        }
        output_stream << "]\n";  // No comma after the last property
    }
    output_stream << "}\n";
}

void write_json_file(const std::string &json_file, Protein &ref, Protein &mem, Dag &d, 
                    std::vector<std::pair<int64_t, int64_t>> &windows, 
                    std::vector<std::pair<int64_t, int64_t>> &windowsp, 
                    std::vector<std::vector<mpq_class>> &ratio,
                    const std::string &alignment_ref,
                    const std::string &alignment_mem)
{
    std::ofstream output_stream;
    output_stream.open(json_file, std::ofstream::app);
    write_json_to_stream(output_stream, ref, mem, d, windows, windowsp, ratio, alignment_ref, alignment_mem);
    output_stream.close();
}

// Overload for backward compatibility
void write_json_file(const std::string &json_file, Protein &ref, Protein &mem, Dag &d, 
                    std::vector<std::pair<int64_t, int64_t>> &windows, 
                    std::vector<std::pair<int64_t, int64_t>> &windowsp, 
                    std::vector<std::vector<mpq_class>> &ratio)
{
    // Call the new version with empty alignment strings
    write_json_file(json_file, ref, mem, d, windows, windowsp, ratio, "", "");
}
