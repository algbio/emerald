#include "cli_handler.h"
#include <getopt.h>
#include <iostream>

void CliHandler::printUsage(const char* programName) {
    const char* usage = R"(
EMERALD - Efficient Method for Evaluating Reliable Alignments
         in Large-scale Datasets

Usage: %s [options]

Required options:
  -f, --input FILE       Input FASTA file with protein sequences
  -o, --output FILE      Output file for safety intervals

Optional parameters:
  -a, --alpha FLOAT      Safety threshold (default: 0.75)
  -d, --delta INT        Suboptimal paths neighborhood (default: 0)
  -g, --gapcost INT      Gap extension cost (default: -1)
  -e, --startgap INT     Gap start cost (default: -11)
  -s, --special INT      Special character score (default: -1)

Output options:
  -j, --json FILE        Write detailed results to JSON file
  -r, --reference STR    Process only alignments with this reference sequence
  -w, --drawgraph DIR    Enable graph visualization, output to directory

Other options:
  -v, --verbose          Enable verbose output
  -h, --help             Show this help message

Example:
  %s -f sequences.fasta -o results.txt -a 0.8
)";
    fprintf(stderr, usage, programName, programName);
}

Config CliHandler::parseCommandLine(int argc, char** argv) {
    Config config;
    
    if (argc < 2) {
        printUsage(argv[0]);
        throw std::runtime_error("No arguments provided. Input and output files are required.");
    }

    static struct option long_options[] = {
        {"alpha", required_argument, nullptr, 'a'},
        {"delta", required_argument, nullptr, 'd'},
        {"input", required_argument, nullptr, 'f'},
        {"output", required_argument, nullptr, 'o'},
        {"costmat", required_argument, nullptr, 'c'},
        {"gapcost", required_argument, nullptr, 'g'},
        {"startgap", required_argument, nullptr, 'e'},
        {"special", required_argument, nullptr, 's'},
        {"json", required_argument, nullptr, 'j'},
        {"reference", required_argument, nullptr, 'r'},
        {"drawgraph", required_argument, nullptr, 'w'},
        {"verbose", no_argument, nullptr, 'v'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    bool has_input = false;
    bool has_output = false;
    while ((opt = getopt_long(argc, argv, "a:d:f:o:c:g:e:s:j:r:w:vh", 
           long_options, nullptr)) != -1) {
        switch (opt) {
            case 'a':
                config.alpha = std::stof(optarg);
                break;
            case 'd':
                config.delta = std::stoul(optarg);
                break;
            case 'f':
                config.input_file = optarg;
                has_input = true;
                break;
            case 'o':
                config.output_file = optarg;
                has_output = true;
                break;
            case 'c':
                config.cost_matrix_file = optarg;
                break;
            case 'g':
                config.gap_cost = std::stoi(optarg);
                break;
            case 'e':
                config.start_gap = std::stoi(optarg);
                break;
            case 's':
                config.special_score = std::stoi(optarg);
                break;
            case 'j':
                config.json_file = optarg;
                break;
            case 'r':
                config.reference_protein = optarg;
                break;
            case 'w':
                config.draw_graph = true;
                config.draw_graph_dir = optarg;
                break;
            case 'v':
                config.verbose = true;
                break;
            case 'h':
                printUsage(argv[0]);
                exit(0);
            default:
                printUsage(argv[0]);
                exit(1);
        }
    }

    if (!has_input || !has_output) {
        printUsage(argv[0]);
        throw std::runtime_error("Missing required option(s): " + 
            std::string(!has_input ? "-f/--input " : "") +
            std::string(!has_output ? "-o/--output" : ""));
    }

    config.validateConfig();
    return config;
}