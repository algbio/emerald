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

		output_stream << "\t\"alignment_graph\": [\n";
		{
			for (int64_t i = 0; i < (int64_t) d.adj.size(); i++) {
				assert (d.adj[i].size() == ratio[i].size());
				output_stream << "\t\t\"" << i << "\": [ ";
				for (int64_t j = 0; j < (int64_t) d.adj[i].size(); j++) {
					output_stream << "\"(" << d.adj[i][j] << "," << ratio[i][j].get_d() << ")\", ";
				}
				output_stream << " ],\n";
			}
		}
		output_stream << "\t],\n";

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
