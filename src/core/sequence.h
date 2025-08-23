#pragma once
#include "sequence_alphabet.h"
#include <string>
#include <memory>

class Sequence {
public:
    Sequence(const std::string& descriptor, const std::string& sequence, 
             std::shared_ptr<SequenceAlphabet> alphabet)
        : descriptor_(descriptor), sequence_(sequence), alphabet_(alphabet) {}
    
    const std::string& getDescriptor() const { return descriptor_; }
    const std::string& getSequence() const { return sequence_; }
    std::shared_ptr<SequenceAlphabet> getAlphabet() const { return alphabet_; }
    
    bool isValid() const {
        for (char c : sequence_) {
            if (!alphabet_->isValidCharacter(c)) return false;
        }
        return true;
    }
    
private:
    std::string descriptor_;
    std::string sequence_;
    std::shared_ptr<SequenceAlphabet> alphabet_;
};