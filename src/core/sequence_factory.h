#pragma once
#include "sequence.h"
#include <memory>
#include <algorithm>

class SequenceFactory {
public:
    static std::unique_ptr<Sequence> createProteinSequence(const std::string& desc, const std::string& seq) {
        return std::make_unique<Sequence>(desc, seq, std::make_shared<ProteinAlphabet>());
    }
    
    static std::unique_ptr<Sequence> createDNASequence(const std::string& desc, const std::string& seq) {
        return std::make_unique<Sequence>(desc, seq, std::make_shared<DNAAlphabet>());
    }
    
    static std::unique_ptr<Sequence> autoDetectSequence(const std::string& desc, const std::string& seq) {
        // Auto-detect based on character composition
        int dna_chars = 0;
        int protein_chars = 0;
        
        for (char c : seq) {
            char upper = std::toupper(c);
            if (upper == 'A' || upper == 'C' || upper == 'G' || upper == 'T') {
                dna_chars++;
            }
            if (upper == 'E' || upper == 'F' || upper == 'I' || upper == 'L' || 
                upper == 'P' || upper == 'Q' || upper == 'Z') {
                protein_chars++; // These are unique to proteins
            }
        }
        
        // If we have protein-specific characters, it's likely protein
        if (protein_chars > 0) {
            return createProteinSequence(desc, seq);
        }
        
        // If >95% ACGT, assume DNA
        if ((double)dna_chars / seq.size() > 0.95) {
            return createDNASequence(desc, seq);
        }
        
        // Default to protein
        return createProteinSequence(desc, seq);
    }
};