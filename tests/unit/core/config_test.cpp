#include <gtest/gtest.h>
#include "../../../src/core/config.h"
#include <filesystem>
#include <fstream>

class ConfigTest : public ::testing::Test {
protected:
    void SetUp() override {
        namespace fs = std::filesystem;
        
        // Create test directories
        fs::create_directories("tests/data/input");
        fs::create_directories("tests/data/output");
        
        // Create a dummy input file
        std::ofstream test_file("tests/data/input/test.fasta");
        test_file.close();
    }

    void TearDown() override {
        std::filesystem::remove_all("tests/data/output");
        std::filesystem::remove_all("tests/data/input");
    }

    Config config;
};

TEST_F(ConfigTest, DefaultsAreInValidRanges) {
    // Test ranges and constraints, not specific values
    EXPECT_GT(config.alpha, 0.5f);
    EXPECT_LE(config.alpha, 1.0f);
    EXPECT_GE(config.delta, 0);
    EXPECT_LT(config.gap_cost, 0);
    EXPECT_LT(config.start_gap, 0);
    EXPECT_LT(config.special_score, 0);
}

TEST_F(ConfigTest, DefaultFlagsAreFalse) {
    EXPECT_FALSE(config.ignore_special);
    EXPECT_FALSE(config.verbose);
    EXPECT_FALSE(config.window_merge);
    EXPECT_FALSE(config.draw_graph);
}

TEST_F(ConfigTest, DefaultPathsAreEmpty) {
    EXPECT_TRUE(config.input_file.empty());
    EXPECT_TRUE(config.output_file.empty());
    EXPECT_TRUE(config.json_file.empty());
    EXPECT_TRUE(config.cost_matrix_file.empty());
    EXPECT_TRUE(config.reference_protein.empty());
}

TEST_F(ConfigTest, CanSetValues) {
    config.alpha = 0.8f;
    config.input_file = "test.fasta";
    config.output_file = "output.txt";
    config.verbose = true;

    EXPECT_FLOAT_EQ(config.alpha, 0.8f);
    EXPECT_EQ(config.input_file, "test.fasta");
    EXPECT_EQ(config.output_file, "output.txt");
    EXPECT_TRUE(config.verbose);
}

TEST_F(ConfigTest, CanSetOptionalPaths) {
    config.json_file = "output.json";
    config.cost_matrix_file = "matrix.txt";
    config.draw_graph_dir = "/tmp/graphs";

    EXPECT_EQ(config.json_file, "output.json");
    EXPECT_EQ(config.cost_matrix_file, "matrix.txt");
    EXPECT_EQ(config.draw_graph_dir, "/tmp/graphs");
}

TEST_F(ConfigTest, ValidatesParameters) {
    // Setup valid base configuration
    config.input_file = "tests/data/input/test.fasta";
    config.output_file = "tests/data/output/result.txt";
    config.alpha = 0.75f;  // valid default
    config.gap_cost = -1;  // valid default
    config.start_gap = -11;  // valid default
    config.special_score = -1;  // valid default
    
    // Test alpha range
    config.alpha = 0.5f;  // boundary
    EXPECT_THROW(config.validateConfig(), std::runtime_error);
    config.alpha = 1.1f;  // too high
    EXPECT_THROW(config.validateConfig(), std::runtime_error);
    config.alpha = 0.75f;  // restore valid value
    EXPECT_NO_THROW(config.validateConfig());

    // Test gap costs
    config.gap_cost = 0;  // invalid
    EXPECT_THROW(config.validateConfig(), std::runtime_error);
    config.gap_cost = -1;  // restore valid value
    EXPECT_NO_THROW(config.validateConfig());

    // Test start gap costs
    config.start_gap = 0;  // invalid
    EXPECT_THROW(config.validateConfig(), std::runtime_error);
    config.start_gap = -11;  // restore valid value
    EXPECT_NO_THROW(config.validateConfig());

    // Test special character score
    config.special_score = 1;  // invalid positive
    EXPECT_THROW(config.validateConfig(), std::runtime_error);
    config.special_score = -1;  // restore valid value
    EXPECT_NO_THROW(config.validateConfig());
}

TEST_F(ConfigTest, ValidatesOptionalFiles) {
    namespace fs = std::filesystem;

    // Setup valid base configuration
    config.input_file = "tests/data/input/test.fasta";
    config.output_file = "tests/data/output/result.txt";

    // Empty optional files should pass both config and file validation
    EXPECT_NO_THROW(config.validateConfig());
    EXPECT_NO_THROW(config.validateFiles());

    // Non-existent cost matrix file should throw from validateFiles()
    config.cost_matrix_file = "nonexistent.mat";
    EXPECT_THROW(config.validateFiles(), std::runtime_error);
    config.cost_matrix_file.clear();  // restore valid state

    // Graph dir: validateFiles does not throw; directories are created by createOutputDirectories()
    config.draw_graph = false;
    config.draw_graph_dir = "tests/data/output/graphs-missing";
    EXPECT_NO_THROW(config.validateFiles());
    EXPECT_NO_THROW(config.createOutputDirectories());  // should be a no-op since draw_graph=false
    EXPECT_FALSE(fs::exists(config.draw_graph_dir));    // not created when draw_graph is false

    // When draw_graph is true, createOutputDirectories() creates the directory
    config.draw_graph = true;
    EXPECT_NO_THROW(config.validateFiles());
    EXPECT_NO_THROW(config.createOutputDirectories());
    EXPECT_TRUE(fs::exists(config.draw_graph_dir));
}