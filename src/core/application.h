#pragma once
#include "config.h"
#include "sequence.h"
#include "../io/sequence_reader.h"
#include "../analysis/alignment_processor.h"
#include <vector>
#include <memory>

class Application {
public:
    explicit Application(
        const Config& config,
        std::unique_ptr<SequenceReader> reader = std::make_unique<DefaultSequenceReader>()
    ) : config_(config),
        alignment_processor_(config),
        sequence_reader_(std::move(reader)) {}

    int run() {
        try {
            auto sequences = readSequences();
            alignment_processor_.processSequences(sequences);
            return 0;
        } catch (const std::exception& e) {
            std::cerr << "Application error: " << e.what() << "\n";
            return 1;
        }
    }

private:
    std::vector<std::unique_ptr<Sequence>> readSequences() {
        auto sequences = sequence_reader_->readFasta(
            config_.input_file, 
            config_.sequence_type
        );
        
        if (sequences.empty()) {
            throw std::runtime_error("No valid sequences found in input file");
        }

        return sequences;
    }

    const Config& config_;
    AlignmentProcessor alignment_processor_;
    std::unique_ptr<SequenceReader> sequence_reader_;
};