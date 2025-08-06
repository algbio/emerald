#pragma once
#include "../core/sequence.h"
#include <vector>
#include <memory>

class SequenceReader {
public:
    static std::vector<std::unique_ptr<Sequence>> readSequences(
        const std::string& filename,
        bool ignoreSpecial
    ) {
        // Stub implementation for testing
        return {};
    }
};