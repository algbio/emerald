#pragma once
#include <string>
#include <vector>

class Sequence {
protected:
    std::string sequence_;
    std::string descriptor_;

public:
    Sequence(const std::string& descriptor, const std::string& sequence)
        : sequence_(sequence), descriptor_(descriptor) {}
    
    virtual ~Sequence() = default;
    
    const std::string& getSequence() const { return sequence_; }
    const std::string& getDescriptor() const { return descriptor_; }
    
    virtual bool isValidSymbol(char symbol) const = 0;
    virtual std::vector<char> getAlphabet() const = 0;
};