#pragma once
#include <string>
#include <unordered_map>
#include <memory>

enum class SequenceType {
    PROTEIN,
    DNA,
    RNA,
    CUSTOM
};

class SequenceAlphabet {
public:
    virtual ~SequenceAlphabet() = default;
    virtual bool isValidCharacter(char c) const = 0;
    virtual size_t getCharacterIndex(char c) const = 0;
    virtual size_t getAlphabetSize() const = 0;
    virtual SequenceType getType() const = 0;
    virtual std::string getName() const = 0;
};

class ProteinAlphabet : public SequenceAlphabet {
public:
    bool isValidCharacter(char c) const override;
    size_t getCharacterIndex(char c) const override;
    size_t getAlphabetSize() const override { return 20; }
    SequenceType getType() const override { return SequenceType::PROTEIN; }
    std::string getName() const override { return "Protein"; }
    
private:
    static const std::unordered_map<char, size_t> amino_acid_table_;
};

class DNAAlphabet : public SequenceAlphabet {
public:
    bool isValidCharacter(char c) const override;
    size_t getCharacterIndex(char c) const override;
    size_t getAlphabetSize() const override { return 4; }
    SequenceType getType() const override { return SequenceType::DNA; }
    std::string getName() const override { return "DNA"; }
    
private:
    static const std::unordered_map<char, size_t> nucleotide_table_;
};