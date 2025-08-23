#pragma once
#include "../core/sequence_alphabet.h"
#include <vector>
#include <string>
#include <memory>

class ScoringMatrix {
public:
    // Default constructor - use BLOSUM62 for proteins
    ScoringMatrix() 
        : alphabet_(std::make_shared<ProteinAlphabet>()),
          special_score_(-5) {
        initializeDefaultMatrix();
    }
    
    // Constructor with file
    ScoringMatrix(const std::string& filename) 
        : alphabet_(std::make_shared<ProteinAlphabet>()),
          special_score_(-5) {
        loadFromFile(filename);
    }
    
    // Constructor with specific alphabet
    ScoringMatrix(std::shared_ptr<SequenceAlphabet> alphabet) 
        : alphabet_(alphabet),
          special_score_(-5) {
        initializeDefaultMatrix();
    }
    
    int64_t getScore(char a, char b) const {
        if (!isValidCharacter(a) || !isValidCharacter(b)) {
            return special_score_;
        }
        
        size_t idx_a = alphabet_->getCharacterIndex(a);
        size_t idx_b = alphabet_->getCharacterIndex(b);
        return matrix_[idx_a][idx_b];
    }
    
    bool isValidCharacter(char c) const {
        return alphabet_->isValidCharacter(c);
    }
    
    void setSpecialScore(int64_t score) {
        special_score_ = score;
    }
    
    SequenceType getAlphabetType() const {
        return alphabet_->getType();
    }
    
    std::shared_ptr<SequenceAlphabet> getAlphabet() const {
        return alphabet_;
    }
    
private:
    std::shared_ptr<SequenceAlphabet> alphabet_;
    std::vector<std::vector<int64_t>> matrix_;
    int64_t special_score_;
    
    void initializeDefaultMatrix();
    void loadFromFile(const std::string& filename);
};