#include "scoring_matrix.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <assert.h>

void ScoringMatrix::initializeDefaultMatrix() {
    size_t alphabet_size = alphabet_->getAlphabetSize();
    matrix_.assign(alphabet_size, std::vector<int64_t>(alphabet_size, 0));
    
    if (alphabet_->getType() == SequenceType::PROTEIN) {
        // BLOSUM62 matrix
        int64_t blosum62[20][20] = {
            // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
            {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0 }, // A
            { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3 }, // R
            { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3 }, // N
            { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3 }, // D
            {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1 }, // C
            { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2 }, // Q
            { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2 }, // E
            {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3 }, // G
            { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3 }, // H
            { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3 }, // I
            { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1 }, // L
            { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2 }, // K
            { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1 }, // M
            { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1 }, // F
            { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2 }, // P
            {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2 }, // S
            {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0 }, // T
            { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3 }, // W
            { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1 }, // Y
            {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4 }  // V
        };
        
        for (size_t i = 0; i < alphabet_size; ++i) {
            for (size_t j = 0; j < alphabet_size; ++j) {
                matrix_[i][j] = blosum62[i][j];
            }
        }
        for (size_t i = 0; i < alphabet_size; ++i) {
            for (size_t j = i + 1; j < alphabet_size; ++j) {
                assert (matrix_[j][i] == matrix_[i][j]); // Ensure symmetry
            }
        }
    } else if (alphabet_->getType() == SequenceType::DNA) {
        // Simple identity matrix for DNA
        for (size_t i = 0; i < alphabet_size; ++i) {
            for (size_t j = 0; j < alphabet_size; ++j) {
                matrix_[i][j] = (i == j) ? 0 : -1;  // edit distance: 0 for match, -1 for mismatch
            }
        }
    }
}

void ScoringMatrix::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open scoring matrix file: " + filename);
    }
    
    size_t alphabet_size = alphabet_->getAlphabetSize();
    matrix_.assign(alphabet_size, std::vector<int64_t>(alphabet_size, 0));
    
    std::string line;
    int row = 0;
    
    while (std::getline(file, line) && row < alphabet_size) {
        std::istringstream iss(line);
        int col = 0;
        int64_t value;
        
        while (iss >> value && col < alphabet_size) {
            matrix_[row][col] = matrix_[col][row] = value;
            ++col;
        }
        
        if (col != row + 1) {
            throw std::runtime_error("Invalid scoring matrix format at row " + std::to_string(row));
        }
        ++row;
    }
    
    if (row != alphabet_size) {
        throw std::runtime_error("Scoring matrix must be " + std::to_string(alphabet_size) + "x" + std::to_string(alphabet_size));
    }
}