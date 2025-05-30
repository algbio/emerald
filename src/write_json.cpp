#include <vector>
#include <string>
#include <fstream>
#include <gmpxx.h>

#include "write_json.h"

void write_json_file(const std::string &json_file, Protein &ref, Protein &mem, Dag &d, std::vector<std::pair<int64_t, int64_t>> &windows, std::vector<std::pair<int64_t, int64_t>> &windowsp, std::vector<std::vector<mpq_class>> &ratio) {
    std::ofstream output_stream;
    output_stream.open(json_file, std::ofstream::app);

    output_stream << "{\n";
    {
        output_stream << "\t\"representative_descriptor\": \"" << ref.descriptor << "\",\n";
        output_stream << "\t\"representative_string\": \"" << ref.sequence << "\",\n";
        output_stream << "\t\"mem_descriptor\": \"" << mem.descriptor << "\",\n";
        output_stream << "\t\"member_string\": \"" << mem.sequence << "\",\n";

        output_stream << "\t\"alignment_graph\": {\n";
        {
            bool first = true;
            for (int64_t i = 0; i < (int64_t) d.adj.size(); i++) {
                assert (d.adj[i].size() == ratio[i].size());
                
                if (!first) output_stream << ",\n";
                first = false;
                
                // Get the index pair for this node
                auto coords = d.transr.at(i);
                int64_t ref_pos = coords.first;
                int64_t mem_pos = coords.second;
                
                output_stream << "\t\t\"(" << ref_pos << "," << mem_pos << ")\": [ ";
                for (int64_t j = 0; j < (int64_t) d.adj[i].size(); j++) {
                    // Get the index pair for the adjacent node
                    auto adj_coords = d.transr.at(d.adj[i][j]);
                    int64_t adj_ref_pos = adj_coords.first;
                    int64_t adj_mem_pos = adj_coords.second;
                    
                    output_stream << "\"((" << adj_ref_pos << "," << adj_mem_pos << ")," << ratio[i][j].get_d() << ")\", ";
                }
                output_stream << " ]";
            }
            output_stream << "\n\t},\n";
        }

        output_stream << "\t\"windows_representative\": [ ";
            for (auto [l, r]: windows) output_stream << "\"(" << l << "," << r << ")\", ";
        output_stream << "],\n";
        output_stream << "\t\"windows_member\": [ ";
            for (auto [l, r]: windowsp) output_stream << "\"(" << l << "," << r << ")\", ";
        output_stream << "],\n";

    }
    output_stream << "}\n";

    output_stream.close();
}
