#include "sequence_alphabet.h"
#include <cctype>

// Initialize protein alphabet mapping
const std::unordered_map<char, size_t> ProteinAlphabet::amino_acid_table_ = {
    {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4}, 
    {'Q', 5}, {'E', 6}, {'G', 7}, {'H', 8}, {'I', 9}, 
    {'L', 10}, {'K', 11}, {'M', 12}, {'F', 13}, {'P', 14}, 
    {'S', 15}, {'T', 16}, {'W', 17}, {'Y', 18}, {'V', 19}
};

// Initialize DNA alphabet mapping
const std::unordered_map<char, size_t> DNAAlphabet::nucleotide_table_ = {
    {'A', 0}, {'C', 1}, {'G', 2}, {'T', 3}
};

bool ProteinAlphabet::isValidCharacter(char c) const {
    char upper_c = std::toupper(c);
    return amino_acid_table_.find(upper_c) != amino_acid_table_.end();
}

size_t ProteinAlphabet::getCharacterIndex(char c) const {
    char upper_c = std::toupper(c);
    auto it = amino_acid_table_.find(upper_c);
    if (it != amino_acid_table_.end()) {
        return it->second;
    }
    throw std::runtime_error("Invalid amino acid: " + std::string(1, c));
}

bool DNAAlphabet::isValidCharacter(char c) const {
    char upper_c = std::toupper(c);
    return nucleotide_table_.find(upper_c) != nucleotide_table_.end();
}

size_t DNAAlphabet::getCharacterIndex(char c) const {
    char upper_c = std::toupper(c);
    auto it = nucleotide_table_.find(upper_c);
    if (it != nucleotide_table_.end()) {
        return it->second;
    }
    throw std::runtime_error("Invalid nucleotide: " + std::string(1, c));
}