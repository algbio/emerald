#pragma once
#include "../core/config.h"
#include "../core/sequence.h"
#include "../io/alignment_writer.h"
#include <vector>
#include <memory>
#include <random>
#include <numeric>

class AlignmentProcessor {
public:
    explicit AlignmentProcessor(const Config& config)
        : config_(config),
          alignment_writer_(config) {}

    // Main processing method
    void processSequences(const std::vector<std::unique_ptr<Sequence>>& sequences) {
        if (sequences.empty()) {
            throw std::runtime_error("No sequences to process");
        }

        if (!config_.reference_protein.empty()) {
            processWithReference(sequences);
        } else {
            processAllPairs(sequences);
        }
    }

private:
    // Core computation methods
    void processWithReference(const std::vector<std::unique_ptr<Sequence>>& sequences) {
        size_t ref_idx = findReferenceSequence(sequences);
        
        for (size_t i = 0; i < sequences.size(); ++i) {
            if (i != ref_idx) {
                computeAlignment(*sequences[ref_idx], *sequences[i]);
            }
        }
    }

    void processAllPairs(const std::vector<std::unique_ptr<Sequence>>& sequences) {
        std::vector<size_t> indices(sequences.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        if (!config_.verbose) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::shuffle(indices.begin(), indices.end(), gen);
        }

        for (size_t i = 0; i < sequences.size(); ++i) {
            for (size_t j = i + 1; j < sequences.size(); ++j) {
                computeAlignment(*sequences[indices[i]], *sequences[indices[j]]);
            }
        }
    }

    size_t findReferenceSequence(const std::vector<std::unique_ptr<Sequence>>& sequences) {
        for (size_t i = 0; i < sequences.size(); ++i) {
            if (sequences[i]->getDescriptor() == config_.reference_protein) {
                return i;
            }
        }
        throw std::runtime_error("Reference protein not found: " + config_.reference_protein);
    }

    void computeAlignment(const Sequence& seq1, const Sequence& seq2);  // Implementation in .cpp

    const Config& config_;
    AlignmentWriter alignment_writer_;
};