#include "alignment_writer.h"
#include <fstream>
#include <iostream>
#include <filesystem>

AlignmentWriter::AlignmentWriter(const Config& config) : config_(config) {}

void AlignmentWriter::writeAlignmentResults(
    const std::vector<std::pair<size_t, size_t>>& windows,
    const Sequence& ref,
    const Sequence& query
) {
    std::ofstream out(config_.output_file);
    if (!out) {
        throw std::runtime_error("Could not open output file: " + config_.output_file);
    }

    out << "# Alignment Results\n";
    out << "# Reference: " << ref.getDescriptor() << "\n";
    out << "# Query: " << query.getDescriptor() << "\n";
    
    for (const auto& window : windows) {
        out << window.first << "\t" << window.second << "\n";
    }
}

void AlignmentWriter::writeJsonReport(
    const std::vector<std::pair<size_t, size_t>>& windows,
    const Sequence& ref,
    const Sequence& query
) {
    if (config_.json_file.empty()) {
        return;
    }

    std::ofstream json(config_.json_file);
    if (!json) {
        throw std::runtime_error("Could not open JSON file: " + config_.json_file);
    }

    json << "{\n"
         << "  \"reference\": \"" << ref.getDescriptor() << "\",\n"
         << "  \"query\": \"" << query.getDescriptor() << "\",\n"
         << "  \"windows\": [\n";
    
    for (size_t i = 0; i < windows.size(); ++i) {
        json << "    {\"start\": " << windows[i].first 
             << ", \"end\": " << windows[i].second << "}";
        if (i < windows.size() - 1) json << ",";
        json << "\n";
    }
    
    json << "  ]\n}\n";
}

void AlignmentWriter::writeAlignmentGraph(
    const std::vector<std::pair<size_t, size_t>>& paths,
    const Sequence& ref,
    const Sequence& query
) {
    if (!config_.draw_graph) {
        return;
    }

    namespace fs = std::filesystem;
    fs::create_directories(config_.draw_graph_dir);
    
    std::string filename = config_.draw_graph_dir + "/" + 
                          ref.getDescriptor() + "_" + 
                          query.getDescriptor() + ".dot";
    
    std::ofstream dot(filename);
    if (!dot) {
        throw std::runtime_error("Could not open graph file: " + filename);
    }

    dot << "digraph alignment {\n"
        << "  rankdir=LR;\n"
        << "  node [shape=circle];\n";
    
    for (const auto& path : paths) {
        dot << "  " << path.first << " -> " << path.second << ";\n";
    }
    
    dot << "}\n";
}