#pragma once
#include "../core/sequence.h"
#include "../core/sequence_factory.h"
#include <string>
#include <vector>
#include <memory>

class SequenceReader {
public:
    virtual ~SequenceReader() = default;
    
    virtual std::vector<std::unique_ptr<Sequence>> readFasta(
        const std::string& filename, 
        SequenceType type = SequenceType::PROTEIN) = 0;
};

class DefaultSequenceReader : public SequenceReader {
public:
    std::vector<std::unique_ptr<Sequence>> readFasta(
        const std::string& filename, 
        SequenceType type = SequenceType::PROTEIN) override;
};