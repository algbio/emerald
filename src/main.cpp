//#include <cassert>
#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <string>
#include <sstream>
#include <sys/stat.h>

#include <gmpxx.h>
#include <omp.h>

#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"
#include "draw_subgraph.h"

int64_t print_usage(char **argv, int64_t help) {
	std::cout << "Usage: " << argv[0] << " -f <clusterfile> -o <outputfile> [arguments]\n\n";
	std::cout << "Optional arguments:\n";
	std::cout << "\t-a, --alpha <value>      \tFloating value, choose edges that appear in (alpha*100)% of all\n";
	std::cout << "\t                         \t(sub-)optimal paths to be safe. (Default: 0.75, Range: (0.5, 1])\n";
	std::cout << "\t-d, --delta <value>      \tInteger value, defines suboptimal paths to be in the delta neighborhood of\n";
	std::cout << "\t                         \tthe optimal. (Default: 0, Range: [0.0, inf))\n";
	std::cout << "\t-c, --costmat <file>     \tReads the aligning score of two symbols from a text file.\n";
	std::cout << "\t                         \tThe text file is a lower triangular matrix with 20 lines. (Default: BLOSUM62)\n";
	std::cout << "\t-g, --gapcost <value>    \tInteger, set the score of aligning a character to a gap. (Default: 1)\n";
	std::cout << "\t-e, --startgap <value>   \tInteger, set the score of starting a gap alignment. (Default: 11)\n";
	std::cout << "\t-s, --special <value>    \tInteger, sets the score of aligning symbols with special characters.\n";
	std::cout << "\t                         \tINF value ignores these charachters. (Default: 1)\n";
	std::cout << "\t-i, --threads <value>    \tInteger, specifies the number of threads (Default: 1).\n";
	std::cout << "\t-r, --reference <protein>\tProtein identity, selects reference protein. By default, there is no reference, and every pair of\n";
	std::cout << "\t                         \tproteins is aligned with each other.\n";
	std::cout << "\t-w, --drawgraph <dir>    \tReturns dot code files of all alignments in an existent directory for plotting\n";
	std::cout << "\t                         \tthe Delta suboptimal subgraph.\n";
	std::cout << "\t-k, --alignments <file> \tCreates a fasta file in the current working directory containing one alignment for every\n";
	std::cout << "\t                         \tsequence pair. Each alignment traverses every safety window. (Default: None)\n";
	std::cout << "\t-m, --windowmerge        \tMerge safety windows if they intersect or are adjacent. EMERALD prints both merged and\n";
	std::cout << "\t                         \tunmerged safety windows if this option is in use. (Default: Off)\n";
	std::cout << "\t-h, --help               \tShows this help message.\n";
	return help;
}

bool file_exists(const std::string &filename) {
	struct stat buffer;
	return (stat (filename.c_str(), &buffer) == 0);
}

static int verbose_flag;
bool drawgraph = false;
std::string drawgraph_dir = ".";
bool window_merge = false;

float alpha = 0.75, TH = 0;
int64_t delta = 0;

int64_t GAP_COST = -1;
int64_t START_GAP = -11;
int64_t SP = -1;
bool ignore_special = false;
std::string input_file, output_file, cost_matrix_file, file_without_path, file_without_path_and_ending, print_alignments;
bool help_flag = false;
bool input_exists = false, output_exists = false;
bool read_cost_matrix = false;
std::string reference = "";
int64_t threads = 1;

// BLOSUM62 matrix
int64_t cost_matrix[21][21] = {
	// Ala  Arg  Asn  Asp  Cys  Gln  Glu  Gly  His  Ile  Leu  Lys  Met  Phe  Pro  Ser  Thr  Trp  Tyr  Val  Def
	{   4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0,  SP }, // Ala
	{  -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3,  SP }, // Arg
	{  -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3,  SP }, // Asn
	{  -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3,  SP }, // Asp
	{   0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1,  SP }, // Cys
	{  -1,   1,   0,   0,  -3,   5,   2,  -2,   0,  -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2,  SP }, // Gln
	{  -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2,  SP }, // Glu
	{   0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3,  SP }, // Gly
	{  -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3,  SP }, // His
	{  -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,  -0,  -3,  -2,  -1,  -3,  -1,   3,  SP }, // Ile
	{  -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4,  -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1,  SP }, // Leu
	{  -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,  SP }, // Lys
	{  -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1,  SP }, // Met
	{  -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1,  SP }, // Phe
	{  -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2,  SP }, // Pro
	{   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2,  SP }, // Ser
	{   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0,  SP }, // Thr
	{  -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2,  11,   2,  -3,  SP }, // Trp
	{  -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1,  SP }, // Tyr
	{   0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4,  SP }, // Val
	{  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP,  SP }, // Def
};

std::vector<Protein> proteins;
int64_t global_ref = -1; // reference protein
std::vector<int64_t> random_order_of_alignments;

void run_case(const int64_t j, const int64_t ref, std::vector<std::stringbuf> &output) {
	int64_t i = random_order_of_alignments[j];
	std::ostream output_stream(&(output[i]));
	const std::string &a = proteins[ref].sequence;
	const std::string &b = proteins[i].sequence;
	output_stream << proteins[i].descriptor << '\n' << b << '\n';

	// Create suboptimal space
	bool random_alignment_as_optimal = false; // TODO: command line option
	Dag d = gen_dag(a, b, cost_matrix, delta, GAP_COST, START_GAP, random_alignment_as_optimal, verbose_flag);
	std::vector<std::vector<int64_t>> adj = d.adj;

	// Find the number of v-t paths (am) and the number of s-v paths (ram) for all nodes v
	std::vector<std::vector<int64_t>> radj((int64_t) adj.size());
	for (int64_t i = 0; i < (int64_t) adj.size(); i++)
		for (int64_t v: adj[i]) radj[v].push_back(i);
	std::vector<mpz_class> am = number_of_paths(adj);
	std::vector<mpz_class> ram = number_of_paths(radj);

	// Find a path that contains all safety windows (it always exists if 0.5 < alpha <= 1)
	std::vector<std::vector<mpq_class>> ratios = path_ratios(d, am, ram);
	std::vector<int64_t> path = find_alpha_path(d, ratios, alpha, verbose_flag);

	if (print_alignments != "") {
		alignment_with_safety_windows(d, path, proteins[ref], proteins[i], print_alignments);
	}
	
	// Make sure that every node is contained only once in the path
	std::unordered_map<int64_t, int64_t> cnt;
	for (int64_t v: path) {
		assert(cnt[v] == 0);
		cnt[v]++;
	}

	// Find safety windows
	std::vector<mpq_class> r = find_ratios(path, adj, ratios);
	auto [swindows, window_ratios, number_of_edges] = safety_windows(am, ram, path, alpha);

	if (drawgraph) {
		std::string dot = draw_subgraph(i, (int64_t) a.size() + 1, (int64_t) b.size() + 1, d, ratios, alpha, a, b);
		std::string file_g = drawgraph_dir + "/dotcode_fasta_" + std::to_string(i) + ".dot";
		std::ofstream str(file_g, std::ofstream::out);
		str << dot;
		str.close();
	}

	std::vector<std::pair<int64_t, int64_t>> windows, windowsp;
	
	std::map<int64_t, std::pair<int64_t, int64_t>> transr = d.transr;
	for (int64_t i = 0; i < (int64_t) swindows.size(); i++) {
		auto [LT, RT] = swindows[i];
		int64_t L = transr[LT].first, R = transr[RT].first;
		int64_t Lp = transr[LT].second, Rp = transr[RT].second;

		windows.emplace_back(L, R);
		windowsp.emplace_back(Lp, Rp);
	}

	auto merge_window = [&](std::vector<std::pair<int64_t, int64_t>> &merged_intervals, std::pair<int64_t, int64_t> merged) {
		if (!merged_intervals.empty() && merged_intervals.back().second >= merged.first) {
			merged.first = merged_intervals.back().first;
			merged_intervals.pop_back();
		}
		merged_intervals.push_back(merged);
	};
	std::vector<std::pair<int64_t, int64_t>> merged_intervals_ref, merged_intervals_member;
	output_stream << windows.size() << '\n';
	for (int64_t k = 0; k < (int64_t) windows.size(); k++) {
		auto [x, y] = windows[k];
		auto [xp, yp] = windowsp[k];
		output_stream << x << ' ' << y << ' ' << xp << ' ' << yp << '\n';

		merge_window(merged_intervals_ref, std::make_pair(x, y));
		merge_window(merged_intervals_member, std::make_pair(xp, yp));
	}
	auto print_sequences_with_sw_coloured = [&](const std::string &s, const std::vector<std::pair<int64_t, int64_t>> &merged_intervals) {
		for (int64_t ch = 0, k = 0; ch < (int64_t) s.size(); ch++) {
			int64_t x = merged_intervals[k].first;
			int64_t y = merged_intervals[k].second;
			while (ch >= y) {
				if (k+1 == (int64_t) merged_intervals.size()) break;
				k++;
				x = merged_intervals[k].first;
				y = merged_intervals[k].second;
			}
			if (x <= ch && ch < y) std::cout << "\033[1;42m" << s[ch] << "\033[0m";
			else std::cout << s[ch];
		}
		std::cout << '\n';
	};
	std::cout << "Alignment (" << j+1 << '/' << ref << ")\n" << proteins[ref].descriptor << '\n' << proteins[i].descriptor << '\n';
	print_sequences_with_sw_coloured(a, merged_intervals_ref);
	print_sequences_with_sw_coloured(b, merged_intervals_member);
	std::cout << '\n';

	if (window_merge) {
		output_stream << "Merged representative safety windows: " << merged_intervals_ref.size() << '\n';
		for (auto [x, y]: merged_intervals_ref)
			output_stream << x << ' ' << y << '\n';
		output_stream << "Merged cluster member safety windows: " << merged_intervals_member.size() << '\n';
		for (auto [xp, yp]: merged_intervals_member)
			output_stream << xp << ' ' << yp << '\n';
	}

	if (verbose_flag) {
		int64_t safe_positions = 0;
		for (auto [L, R]: merged_intervals_member) safe_positions += R - L;
		std::vector<std::pair<int64_t, int64_t>> sno = safe_not_opt(path, swindows, d, false);
		std::cout << "NUM OF EDGES: " << number_of_edges << std::endl;
		std::cout << "SAFE POSITIONS OF CLUSTER MEMBER: " << safe_positions << std::endl;
		output_stream << "Safe edges not included in optimal paths: " << sno.size() << " (" << double(sno.size())/double(number_of_edges) * 100 << "%)\n";
	}
}

int main(int argc, char **argv) {
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output


	int64_t c;
	while (1) {
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 },
			{ "alpha", required_argument, 0, 'a' },
			{ "delta", required_argument, 0, 'd' },
			{ "costmat", required_argument, 0, 'c' },
			{ "gapcost", required_argument, 0, 'g' },
			{ "startgap", required_argument, 0, 'e' },
			{ "special", required_argument, 0, 's' },
			{ "file", required_argument, 0, 'f' },
			{ "output", required_argument, 0, 'o' },
			{ "help", no_argument, 0, 'h' },
			{ "threads", required_argument, 0, 'i' },
			{ "reference", required_argument, 0, 'r' },
			{ "drawgraph", required_argument, 0, 'w' },
			{ "alignments", required_argument, 0, 'k' },
			{ "windowmerge", no_argument, 0, 'm' },
			{ 0, 0, 0, 0 }
		};
	
		int option_index = 0;
		c = getopt_long(argc, argv, "a:d:c:g:e:s:f:o:hi:r:wk:m", long_options, &option_index);
		if (c == -1) break;

		switch (c) {
			case 0:
				if (long_options[option_index].flag != 0)
					break;
				std::cout << "Option " << long_options[option_index].name;
				if (optarg)
					std::cout << " with arg " << optarg;
				std::cout << std::endl;
				break;
			case 'a':
				alpha = std::stof(optarg);
				break;
			case 'd':
				delta = std::stof(optarg);
				break;
			case 'c':
				read_cost_matrix = true;
				cost_matrix_file = optarg;
				break;
			case 'g':
				GAP_COST = atoi(optarg);
				break;
			case 'e':
				START_GAP = atoi(optarg);
				break;
			case 's':
				if (strcmp(optarg, "INF") == 0) ignore_special = true;
				else SP = atoi(optarg);
				break;
			case 'f':
				input_exists = true;
				input_file = optarg;
				file_without_path = input_file.substr(input_file.find_last_of("/\\") + 1);
				if (file_without_path == "") file_without_path = input_file; // input_file was in the same directory and does not contain / or '\\'
				file_without_path_and_ending = file_without_path.substr(0, file_without_path.find_last_of("."));
				break;
			case 'o':
				output_exists = true;
				output_file = optarg;
				break;
			case 'h':
				help_flag = true;
				break;
			case 'i':
				threads = atoi(optarg);
				break;
			case 'r':
				reference = optarg;
				break;
			case 'w':
				drawgraph_dir = optarg;
				drawgraph = true;
				break;
			case 'k':
				print_alignments = optarg;
				break;
			case 'm':
				window_merge = true;
			case '?': break;
			default: abort();
		}
	}
	
	if (help_flag)
		return print_usage(argv, 0);

	if (!input_exists || !output_exists) {
		std::cerr << "Error: No " << (!input_exists ? "input" : "output") << " file given.\n";
		return print_usage(argv, 1);
	}

	if (drawgraph && !file_exists(drawgraph_dir)) {
		std::cerr << "Error: directory " << drawgraph_dir << " does not exist.\n";
		return 2;
	}

	// Welcome message
	std::cout << "EMERALD is free software (GPL v3) and is maintained on https://github.com/algbio/emerald.\n";
	std::cout << "Please cite the following reference when using EMERALD for your research:\n";
	std::cout << "Grigorjew, A., Gynter, A., Dias, F.H. et al. Sensitive inference of alignment-safe intervals from biodiverse protein sequence clusters using EMERALD. Genome Biol 24, 168 (2023). https://doi.org/10.1186/s13059-023-03008-6\n\n";

	if (alpha > 1.0)
		std::cerr << "Warning: for alpha values > 1.0 the program will not return any safety windows.\n";
	else if (alpha <= 0.5)
		std::cerr << "Warning: for alpha values <= 0.5 the program will not behave well defined and might crash.\n";

	if (file_exists(print_alignments)) {
		std::cerr << "Warning: alignment file " << print_alignments << " already exists. The alignments will be appended to the file.\n";
	}

	if (read_cost_matrix) {
		std::ifstream costmat(cost_matrix_file);
		for (int64_t i = 0; i < 20; i++) for (int64_t j = 0; j <= i; j++) {
			costmat >> cost_matrix[i][j];
			cost_matrix[j][i] = cost_matrix[i][j];
		}
	}

	// TODO: Read these from the LTA array instead
	std::vector<char> special_chars = { 'B', 'X', 'Z' };
	auto contains_special = [&](const std::string &a) {
		for (char c: special_chars)
			if (a.find(c) != std::string::npos) return true;
		return false;
	};

	std::ifstream input(input_file);
	proteins.clear();
	for (std::string line; std::getline(input, line); ) {
		if ((int64_t) line.size() <= 0) continue;
		if (line[0] == '>') {
			if ((int64_t) proteins.size() > 0 && ignore_special && contains_special(proteins.back().sequence))
				proteins.pop_back();
			proteins.push_back(Protein(line));
		} else {
			proteins.back().sequence += line;
		}
	}
	if ((int64_t) proteins.size() > 0 && ignore_special && contains_special(proteins.back().sequence))
		proteins.pop_back();

	int64_t PS = (int64_t) proteins.size();

	if (PS == 0) {
		std::cerr << "Error: Protein sequence list is empty.\n";
		return 3;
	}

	// find reference protein
	//if (reference != "") {
	//	bool found = false;
	//	for (int64_t i = 0; i < PS; i++) {
	//		if (proteins[i].descriptor.find(reference) != std::string::npos) {
	//			if (found)
	//				std::cerr << "Reference identity found more than once. Using the last found protein as reference.\n";
	//			ref = i;
	//			found = true;
	//		}
	//	}
	//	if (!found) std::cerr << "Reference identity not found, retreat to aligning all pairs of proteins without reference.\n";
	//}

	std::ofstream output_stream;
	output_stream.open(output_file, std::ofstream::app);

	// reference protein and number of proteins in the cluster
	for (int64_t ref = 1; ref < PS; ref++) {
		output_stream << proteins[ref].descriptor << '\n' << proteins[ref].sequence << '\n';
		std::vector<std::stringbuf> output(PS); // TODO: Do we really want to do this?
		random_order_of_alignments.clear();
		for (int64_t i = 0; i < ref; i++) random_order_of_alignments.push_back(i);
		if (threads > 1) {
			std::random_device rd;
			std::mt19937 g(rd());
			std::shuffle(random_order_of_alignments.begin(), random_order_of_alignments.end(), g);
		}

		#pragma omp parallel for num_threads(threads)
		for (int64_t j = 0; j < ref; j++)
			run_case(j, ref, output);

		for (int64_t i = 0; i < PS; i++) if (i != ref) output_stream << output[i].str();
		if (ref+1 < PS) output_stream << '\n';
	}
	std::cout << "Safety intervals stored in " << output_file << ".\n";
}
