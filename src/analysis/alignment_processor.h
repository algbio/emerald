#pragma once
#include "../core/config.h"
#include "../core/sequence.h"
#include "scoring_matrix.h"
#include "dp_matrix.h"
#include "dag_builder.h"
#include "alignment_result.h"
#include "../io/alignment_writer.h"
#include <vector>
#include <memory>

class AlignmentProcessor {
public:
    using AlignmentResult = ::AlignmentResult;  // Use the shared struct

    explicit AlignmentProcessor(const Config& config) 
        : config_(config), 
          scoring_matrix_(config.cost_matrix_file.empty() ? 
                         ScoringMatrix() : ScoringMatrix(config.cost_matrix_file)),
          alignment_writer_(config) {
        
        // Set special character score from config
        scoring_matrix_.setSpecialScore(config.special_score);
    }

    void processSequences(const std::vector<std::unique_ptr<Sequence>>& sequences);

private:
    void processWithReference(const std::vector<std::unique_ptr<Sequence>>& sequences);
    void processAllPairs(const std::vector<std::unique_ptr<Sequence>>& sequences);
    
    size_t findReferenceSequence(const std::vector<std::unique_ptr<Sequence>>& sequences);
    void computeAlignment(const Sequence& seq1, const Sequence& seq2);
    AlignmentResult computeOptimalAlignment(const std::string& a, const std::string& b);
    
    void validateSequences(const std::string& seq, const std::string& name) const;

    const Config& config_;
    ScoringMatrix scoring_matrix_;
    AlignmentWriter alignment_writer_;
};