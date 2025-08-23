#pragma once
#include "../core/sequence.h"
#include "../core/config.h"
#include "../analysis/alignment_result.h"

class AlignmentWriter {
public:
    using AlignmentResult = ::AlignmentResult;
    
    explicit AlignmentWriter(const Config& config);
    
    // Main method for our new design
    void writeAlignment(
        const Sequence& seq1,
        const Sequence& seq2,
        const AlignmentResult& result
    );

private:
    const Config& config_;
    
    void writeToFile(const std::string& content, const std::string& suffix = "") const;
    void writeJsonOutput(const Sequence& seq1, const Sequence& seq2, const AlignmentResult& result) const;
    void writeTextOutput(const Sequence& seq1, const Sequence& seq2, const AlignmentResult& result) const;
};