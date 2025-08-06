#include <gtest/gtest.h>
#include "../../../src/core/config.h"
#include <filesystem>

class ConfigTest : public ::testing::Test {
protected:
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