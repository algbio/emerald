#include "sequence_reader.h"
#include <fstream>
#include <stdexcept>

std::vector<std::unique_ptr<Sequence>> DefaultSequenceReader::readFasta(
    const std::string& filename, 
    SequenceType type) {
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    std::vector<std::unique_ptr<Sequence>> sequences;
    std::string line, desc, seq;
    bool reading_sequence = false;

    while (std::getline(file, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {  // Description line
            if (reading_sequence && !seq.empty()) {
                // Save previous sequence
                switch (type) {
                    case SequenceType::PROTEIN:
                        sequences.push_back(SequenceFactory::createProteinSequence(desc, seq));
                        break;
                    case SequenceType::DNA:
                        sequences.push_back(SequenceFactory::createDNASequence(desc, seq));
                        break;
                    default:
                        sequences.push_back(SequenceFactory::autoDetectSequence(desc, seq));
                        break;
                }
                seq.clear();
            }
            desc = line.substr(1);  // Remove '>' character
            reading_sequence = true;
        } else if (reading_sequence) {
            // Append to current sequence
            seq += line;
        }
    }

    // Save final sequence
    if (reading_sequence && !seq.empty()) {
        switch (type) {
            case SequenceType::PROTEIN:
                sequences.push_back(SequenceFactory::createProteinSequence(desc, seq));
                break;
            case SequenceType::DNA:
                sequences.push_back(SequenceFactory::createDNASequence(desc, seq));
                break;
            default:
                sequences.push_back(SequenceFactory::autoDetectSequence(desc, seq));
                break;
        }
    }

    // Check if any sequences were found
    if (sequences.empty()) {
        throw std::runtime_error("No valid sequences found in file: " + filename);
    }

    return sequences;
}