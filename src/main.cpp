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

#include <gmpxx.h>
#include <omp.h>

#include "alpha_safe_paths.h"
#include "safety_windows.h"
#include "optimal_paths.h"
#include "draw_subgraph.h"

int64_t print_usage(char **argv, int64_t help) {
	std::cout << "How to run: " << argv[0] << " -f <clusterfile> [OPTION...]\n\n";
	std::cout << "\t-a, --alpha        \tFloating value, choose edges that appear in (alpha*100)% of all\n";
	std::cout << "\t                   \t(sub-)optimal paths to be safe. (Default: 0.75, Range: (0.5, 1])\n";
	std::cout << "\t-p, --approximation\tBoolean, if true, big integers will be approximated by doubles. Output might be less accurate, but running speed will be increased. (Default: false)\n";
	std::cout << "\t-d, --delta        \tInteger value, defines suboptimal paths to be in the delta neighborhood of the optimal. (Default: 0, Range: [0.0, inf))\n";
	std::cout << "\t-c, --costmat      \tReads the aligning score of two symbols from a text file.\n";
	std::cout << "\t                   \tThe text file is a lower triangular matrix with 20 lines. (Default: BLOSUM62)\n";
	std::cout << "\t-g, --gapcost      \tInteger, set the score of aligning a character to a gap. (Default: 1)\n";
	std::cout << "\t-e, --startgap     \tInteger, set the score of starting a gap alignment. (Default: 11)\n";
	std::cout << "\t-s, --special      \tInteger, sets the score of aligning symbols with special characters.\n";
	std::cout << "\t                   \tINF value ignores these charachters. (Default: 1)\n";
	std::cout << "\t-i, --threads      \tInteger, specifies the number of threads (Default: 1).\n";
	std::cout << "\t-r, --reference    \tProtein identity, selects reference protein. By default, this is the first protein.\n";
	std::cout << "\t-w, --drawgraph    \tReturns dot code files of all alignments for plotting the Delta suboptimal subgraph (for debug purposes).\n";
	std::cout << "\t-k, --alignments   \tNon-negative integer n, create a fasta file containing randomly chosen n suboptimal alignments. (Default: 0)\n";
	std::cout << "\t-m, --windowmerge  \tMerge safety windows if they intersect or are adjacent. EMERALD prints both merged and unmerged safety windows if this option is in use. (Default: Off)\n"
	std::cout << "\t-h, --help         \tShows this help message.\n";
	return help;
}

struct Protein {
	std::string descriptor;
	std::string sequence;
	Protein(std::string descriptor) : descriptor(descriptor) {}
};

static int verbose_flag;
bool use_approx = false;
bool drawgraph = false;
bool window_merge = false;

float alpha = 0.75, TH = 0;
int64_t delta = 0;

int64_t GAP_COST = -1;
int64_t START_GAP = -11;
int64_t SP = -1;
bool ignore_special = false;
std::string file, cost_matrix_file, file_without_path;
bool help_flag = false;
bool read_file = false;
bool read_cost_matrix = false;
int64_t print_alignments = 0;
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
int64_t ref = 0; // reference protein
std::vector<int64_t> random_order_of_alignments;

template<class T, class K>
void run_case(const int64_t j, std::vector<std::stringbuf> &output) {
	int64_t i = random_order_of_alignments[j];
	std::ostream output_stream(&(output[i]));
	const std::string &a = proteins[ref].sequence;
	const std::string &b = proteins[i].sequence;
	output_stream << proteins[i].descriptor << '\n' << b << '\n';

	// Create suboptimal space
	bool random_alignment_as_optimal = false; // TODO: command line option
	Dag d = gen_dag<K>(a, b, cost_matrix, delta, GAP_COST, START_GAP, random_alignment_as_optimal, verbose_flag);
	std::vector<std::vector<int64_t>> adj = d.adj;

	if (print_alignments > 0) {
		std::string subop_fasta_file = "suboptimal_" + file_without_path;
		alignments_into_fasta(print_alignments, d, a, b, subop_fasta_file, std::to_string(i), proteins[i].descriptor);
	}

	// Find the number of v-t paths (am) and the number of s-v paths (ram) for all nodes v
	std::vector<std::vector<int64_t>> radj((int64_t) adj.size());
	for (int64_t i = 0; i < (int64_t) adj.size(); i++)
		for (int64_t v: adj[i]) radj[v].push_back(i);
	std::vector<T> am = number_of_paths<T>(adj);
	std::vector<T> ram = number_of_paths<T>(radj);

	// Find a path that contains all safety windows (it always exists if alpha \in (0.5, 1])
	std::vector<std::vector<K>> ratios = path_ratios<T, K>(d, am, ram);
	std::vector<int64_t> path = find_alpha_path<K>(d, ratios, alpha, verbose_flag);
	
	// Just make sure that every node is contained only once in the path
	std::unordered_map<int64_t, int64_t> cnt;
	for (int64_t v: path) {
		assert(cnt[v] == 0);
		cnt[v]++;
	}

	// Find safety windows
	std::vector<K> r = find_ratios<K>(path, adj, ratios);
	auto [swindows, window_ratios, number_of_edges] = safety_windows<T, K>(am, ram, path, alpha);

	if (drawgraph) {
		std::string dot = draw_subgraph<K>(i, (int64_t) a.size() + 1, (int64_t) b.size() + 1, d, ratios, alpha, a, b);
		std::string file_g = "fasta_" + std::to_string(i) + ".dot";
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
	//test_gen_dag();
	//return 0;
	std::cerr << std::fixed << std::setprecision(10); // debug output
	std::cout << std::fixed << std::setprecision(10); // debug output


	int64_t c;
	while (1) {
		static struct option long_options[] = {
			{ "verbose", no_argument, &verbose_flag, 1 }, // this does nothing for now
			{ "alpha", required_argument, 0, 'a' },
			{ "delta", required_argument, 0, 'd' },
			{ "costmat", required_argument, 0, 'c' },
			{ "gapcost", required_argument, 0, 'g' },
			{ "startgap", required_argument, 0, 'e' },
			{ "special", required_argument, 0, 's' },
			{ "file", required_argument, 0, 'f' },
			{ "help", no_argument, 0, 'h' },
			{ "threads", required_argument, 0, 'i' },
			{ "reference", required_argument, 0, 'r' },
			{ "approximation", no_argument, 0, 'p' },
			{ "drawgraph", no_argument, 0, 'w' },
			{ "alignments", required_argument, 0, 'k' },
			{ "windowmerge", no_argument, 0, 'm' },
			{ 0, 0, 0, 0 }
		};
	
		int option_index = 0;
		c = getopt_long(argc, argv, "a:d:c:g:e:s:f:hi:r:pwk:m", long_options, &option_index);
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
				read_file = true;
				file = optarg;
				file_without_path = file.substr(file.find_last_of("/\\") + 1);
				if (file_without_path == "") file_without_path = file; // file was in the same directory and does not contain / or '\\'
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
			case 'p':
				use_approx = true;
				break;
			case 'w':
				drawgraph = true;
				break;
			case 'k':
				print_alignments = atoi(optarg);
				break;
			case 'm':
				window_merge = true;
			case '?': break;
			default: abort();
		}
	}
	
	if (help_flag)
		return print_usage(argv, 0);

	if (!read_file)
		return print_usage(argv, 1);

	if (alpha > 1.0)
		std::cerr << "Warning: for alpha values > 1.0 the program will not return any safety windows.\n";
	else if (alpha <= 0.5)
		std::cerr << "Warning: for alpha values <= 0.5 the program will not behave well defined and might crash.\n";

	if (print_alignments < 0)
		std::cerr << "Warning: alignments value is negative, will be treated as 0.\n";
	else if (print_alignments > 0) {
		std::string subop_fasta_file = "suboptimal_" + file_without_path;
		std::ofstream ofs;
		ofs.open(subop_fasta_file, std::ofstream::out | std::ofstream::trunc);
		ofs.close();
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

	std::ifstream input(file);
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
		std::cout << "Protein sequence list is empty.\n";
		return 2;
	}

	// find reference protein
	if (reference != "") {
		bool found = false;
		for (int64_t i = 0; i < PS; i++) {
			if (proteins[i].descriptor.find(reference) != std::string::npos) {
				if (found)
					std::cerr << "Representative identity found more than once. Using the last found protein as reference.\n";
				ref = i;
				found = true;
			}
		}
		if (!found) std::cerr << "Reference identity not found, using the first protein as reference.\n";
	}

	// reference protein and number of proteins in the cluster
	std::cout << proteins[ref].descriptor << '\n' << proteins[ref].sequence << '\n';
	std::cout << PS << '\n';
	std::vector<std::stringbuf> output(PS);
	random_order_of_alignments.clear();
	for (int64_t i = 0; i < PS; i++) if (i != ref) random_order_of_alignments.push_back(i);
	if (threads > 1) {
		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(random_order_of_alignments.begin(), random_order_of_alignments.end(), g);
	}

	#pragma omp parallel for num_threads(threads)
	for (int64_t j = 0; j < PS - 1; j++)
		if (!use_approx)
			run_case<mpz_class, mpq_class>(j, output);
		else
			run_case<double, double>(j, output);

	for (int64_t i = 0; i < PS; i++) if (i != ref) std::cout << output[i].str();
}
