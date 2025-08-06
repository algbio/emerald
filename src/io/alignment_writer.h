#pragma once
#include "../core/sequence.h"
#include "../core/config.h"
#include <vector>

class AlignmentWriter {
public:
    explicit AlignmentWriter(const Config& config);
    
    void writeAlignmentResults(
        const std::vector<std::pair<size_t, size_t>>& windows,
        const Sequence& ref,
        const Sequence& query
    );

    void writeJsonReport(
        const std::vector<std::pair<size_t, size_t>>& windows,
        const Sequence& ref,
        const Sequence& query
    );

    void writeAlignmentGraph(
        const std::vector<std::pair<size_t, size_t>>& paths,
        const Sequence& ref,
        const Sequence& query
    );

private:
    const Config& config_;
};