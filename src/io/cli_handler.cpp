#include "cli_handler.h"
#include <getopt.h>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <filesystem>

// Helper template for safer numeric parsing with error checking
template <typename T>
T parse_numeric(const char* arg, const std::string& option_name) {
    if (!arg || *arg == '\0') {
        throw std::runtime_error("Missing value for " + option_name);
    }
    
    try {
        std::string arg_str(arg);
        size_t pos = 0;
        
        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>) {
            T value = std::stof(arg_str, &pos);
            if (pos != arg_str.length()) {
                throw std::runtime_error("Invalid characters in " + option_name + " value: " + arg_str);
            }
            return value;
        } 
        else if constexpr (std::is_same_v<T, int64_t> || std::is_same_v<T, int>) {
            T value = std::stoll(arg_str, &pos);
            if (pos != arg_str.length()) {
                throw std::runtime_error("Invalid characters in " + option_name + " value: " + arg_str);
            }
            return value;
        }
        else if constexpr (std::is_same_v<T, uint64_t> || std::is_same_v<T, unsigned int>) {
            T value = std::stoull(arg_str, &pos);
            if (pos != arg_str.length()) {
                throw std::runtime_error("Invalid characters in " + option_name + " value: " + arg_str);
            }
            return value;
        }
    } catch (const std::invalid_argument&) {
        throw std::runtime_error("Invalid format for " + option_name + ": " + arg);
    } catch (const std::out_of_range&) {
        throw std::runtime_error(option_name + " value out of range: " + arg);
    }
    
    return T(); // Shouldn't reach here
}

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
  -d, --delta INT        Suboptimal paths neighborhood (default: 10)
  -g, --gapcost INT      Gap extension cost (default: -1)
  -e, --startgap INT     Gap start cost (default: -11)
  -s, --special VAL      Special character score (default: -5). Use 'INF' to ignore specials.
  -m, --windowmerge      Merge adjacent/intersecting safety windows (default: off)
  -t, --type STR         Sequence type: 'protein' or 'dna' (default: protein)

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
        {"windowmerge", no_argument, nullptr, 'm'},
        {"type", required_argument, nullptr, 't'},
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
    
    while ((opt = getopt_long(argc, argv, "a:d:f:o:c:g:e:s:j:r:w:vmh", 
           long_options, nullptr)) != -1) {
        try {
            switch (opt) {
                case 'a': {
                    config.alpha = parse_numeric<float>(optarg, "alpha");
                    break;
                }
                case 'd': {
                    int64_t delta = parse_numeric<int64_t>(optarg, "delta");
                    config.delta = delta;
                    break;
                }
                case 'f':
                    if (!optarg || *optarg == '\0') {
                        throw std::runtime_error("Empty input file path provided");
                    }
                    config.input_file = optarg;
                    has_input = true;
                    break;
                case 'o':
                    if (!optarg || *optarg == '\0') {
                        throw std::runtime_error("Empty output file path provided");
                    }
                    config.output_file = optarg;
                    has_output = true;
                    break;
                case 'c':
                    if (!optarg || *optarg == '\0') {
                        throw std::runtime_error("Empty cost matrix file path");
                    }
                    config.cost_matrix_file = optarg;
                    break;
                case 'g':
                    config.gap_cost = parse_numeric<int64_t>(optarg, "gap cost");
                    break;
                case 'e':
                    config.start_gap = parse_numeric<int64_t>(optarg, "start gap");
                    break;
                case 's': {
                    // Support "INF" (case-insensitive) to toggle ignore_special
                    std::string sval = optarg ? std::string(optarg) : "";
                    auto to_lower = [](std::string s) {
                        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch){ return std::tolower(ch); });
                        return s;
                    };
                    if (to_lower(sval) == "inf") {
                        config.ignore_special = true;
                    } else {
                        config.special_score = parse_numeric<int64_t>(optarg, "special score");
                    }
                    break;
                }
                case 'm':
                    config.window_merge = true;
                    break;
                case 't': {
                    if (!optarg || *optarg == '\0') {
                        throw std::runtime_error("Empty sequence type provided");
                    }
                    std::string type_str = optarg;
                    if (type_str == "protein") {
                        config.sequence_type = SequenceType::PROTEIN;
                    } else if (type_str == "dna") {
                        config.sequence_type = SequenceType::DNA;
                    } else {
                        throw std::runtime_error("Invalid sequence type: " + type_str + ". Valid options are 'protein' or 'dna'.");
                    }
                    break;
                }
                case 'j':
                    if (!optarg || *optarg == '\0') {
                        throw std::runtime_error("Empty JSON file path");
                    }
                    config.json_file = optarg;
                    break;
                case 'r':
                    if (!optarg || *optarg == '\0') {
                        throw std::runtime_error("Empty reference protein name");
                    }
                    config.reference_protein = optarg;
                    break;
                case 'w':
                    config.draw_graph = true;
                    if (optarg && *optarg) {
                        config.draw_graph_dir = optarg;
                    }
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
        } catch (const std::exception& e) {
            printUsage(argv[0]);
            throw std::runtime_error(std::string("Error parsing arguments: ") + e.what());
        }
    }

    // Check required options - this is the ONLY place we check for required options
    if (!has_input || !has_output) {
        printUsage(argv[0]);
        throw std::runtime_error("Missing required option(s): " + 
            std::string(!has_input ? "-f/--input " : "") +
            std::string(!has_output ? "-o/--output" : ""));
    }

    // Also check that the values aren't empty
    if (config.input_file.empty()) {
        throw std::runtime_error("Empty input file path provided");
    }

    if (config.output_file.empty()) {
        throw std::runtime_error("Empty output file path provided");
    }

    // Validate business rules and create directories
    try {
        config.validateConfig();
        config.validateFiles();
        config.createOutputDirectories();
    } catch (const std::exception& e) {
        printUsage(argv[0]);
        throw;
    }

    return config;
}