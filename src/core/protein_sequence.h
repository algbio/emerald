#pragma once
#include "sequence.h"

class ProteinSequence : public Sequence {
public:
    ProteinSequence(const std::string& descriptor, const std::string& sequence)
        : Sequence(descriptor, sequence) {}
    
    bool isValidSymbol(char symbol) const override {
        static const std::string validAminoAcids = "ACDEFGHIKLMNPQRSTVWY";
        return validAminoAcids.find(symbol) != std::string::npos;
    }
    
    std::vector<char> getAlphabet() const override {
        return {'A','C','D','E','F','G','H','I','K','L',
                'M','N','P','Q','R','S','T','V','W','Y'};
    }
    
    std::vector<char> getSpecialChars() const {
        return {'B', 'X', 'Z'};
    }
};