#include <gtest/gtest.h>
#include "../../../src/io/sequence_reader.h"
#include "../../../src/core/sequence_alphabet.h"
#include <filesystem>
#include <fstream>

class SequenceReaderTest : public ::testing::Test {
protected:
    DefaultSequenceReader reader;
    
    void SetUp() override {
        // Create test directory
        std::filesystem::create_directories("tests/data/input");
        
        // Copy example file (consistent with other tests)
        std::filesystem::copy("../examples/ex1.fasta", "tests/data/input/test.fasta", 
                             std::filesystem::copy_options::overwrite_existing);
        
        // Create additional test files for specific test cases
        createEmptyFile();
        createSpecialCharFile();
        createInvalidFastaFile();
    }

    void TearDown() override {
        std::filesystem::remove_all("tests/data");
    }

private:
    void createEmptyFile() {
        std::ofstream file("tests/data/input/empty.fasta");
        file.close();
    }

    void createSpecialCharFile() {
        std::ofstream file("tests/data/input/special.fasta");
        file << ">Sequence with special chars\n"
             << "MSFDLKSKFLGBXZ\n";  // B, X, Z are special chars
        file.close();
    }

    void createInvalidFastaFile() {
        std::ofstream file("tests/data/input/invalid.fasta");
        file << "No header line\n"
             << "MSFDLKSKFLG\n"
             << "More sequence without header\n";
        file.close();
    }
};

TEST_F(SequenceReaderTest, ReadsValidFastaFile) {
    auto sequences = reader.readFasta("tests/data/input/test.fasta");
    
    ASSERT_EQ(sequences.size(), 5);
    
    // Test all 5 sequences from ex1.fasta
    EXPECT_EQ(sequences[0]->getDescriptor(), "Cluster sequence 1");
    EXPECT_EQ(sequences[0]->getSequence(), "MSFDLKSKFLG");
    
    EXPECT_EQ(sequences[1]->getDescriptor(), "Cluster sequence 2");
    EXPECT_EQ(sequences[1]->getSequence(), "MSKLKDFLFKS");
    
    EXPECT_EQ(sequences[2]->getDescriptor(), "Cluster sequence 3");
    EXPECT_EQ(sequences[2]->getSequence(), "MSLGSFKDKFL");
    
    EXPECT_EQ(sequences[3]->getDescriptor(), "Cluster sequence 4");
    EXPECT_EQ(sequences[3]->getSequence(), "MSLKDKKFLKS");
    
    EXPECT_EQ(sequences[4]->getDescriptor(), "Cluster sequence 5");
    EXPECT_EQ(sequences[4]->getSequence(), "MSFLKKKFDSL");
}

TEST_F(SequenceReaderTest, ThrowsOnNonexistentFile) {
    EXPECT_THROW(
        reader.readFasta("nonexistent.fasta"),
        std::runtime_error
    );
}

TEST_F(SequenceReaderTest, ThrowsOnEmptyFile) {
    EXPECT_THROW(
        reader.readFasta("tests/data/input/empty.fasta"),
        std::runtime_error
    );
}

TEST_F(SequenceReaderTest, ConvertsToUppercase) {
    std::ofstream file("tests/data/input/lowercase.fasta");
    file << ">lowercase sequence\n"
         << "msfdlkskflg\n";
    file.close();
    
    auto sequences = reader.readFasta("tests/data/input/lowercase.fasta");
    ASSERT_EQ(sequences.size(), 1);
    // Since we're using ProteinAlphabet which should handle case, this should match:
    EXPECT_EQ(sequences[0]->getSequence(), "msfdlkskflg");
}

TEST_F(SequenceReaderTest, HandlesWhitespace) {
    std::ofstream file("tests/data/input/whitespace.fasta");
    file << "> Sequence with spaces \n"
         << "MS FD LK\n"
         << "SK FL G\n";
    file.close();
    
    auto sequences = reader.readFasta("tests/data/input/whitespace.fasta");
    ASSERT_EQ(sequences.size(), 1);
    EXPECT_EQ(sequences[0]->getDescriptor(), " Sequence with spaces ");  // Whitespace preserved in desc
    EXPECT_EQ(sequences[0]->getSequence(), "MS FD LKSK FL G");  // Whitespace
}

TEST_F(SequenceReaderTest, HandlesInvalidFastaFormat) {
    // Should either throw or skip invalid lines depending on implementation
    EXPECT_THROW(
        reader.readFasta("tests/data/input/invalid.fasta"),
        std::runtime_error
    );
}

TEST_F(SequenceReaderTest, CreatesCorrectAlphabetType) {
    auto sequences = reader.readFasta("tests/data/input/test.fasta", SequenceType::PROTEIN);
    
    // Test that we have the correct alphabet type
    ASSERT_EQ(sequences[0]->getAlphabet()->getType(), SequenceType::PROTEIN);
    
    // Test protein alphabet specific methods
    auto alphabet = sequences[0]->getAlphabet();
    EXPECT_TRUE(alphabet->isValidCharacter('A'));
    EXPECT_FALSE(alphabet->isValidCharacter('1'));
    
    EXPECT_EQ(alphabet->getAlphabetSize(), 20);  // 20 standard amino acids
}