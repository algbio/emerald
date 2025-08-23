#include <gtest/gtest.h>
#include "../../../src/io/cli_handler.h"
#include "../../../src/core/config.h"
#include <string>
#include <cstdio>
#include <filesystem>
#include <getopt.h> 

class CliHandlerTest : public ::testing::Test {
protected:
    void SetUp() override {
        namespace fs = std::filesystem;
        
        // Reset getopt state for each test
        optind = 1;
        opterr = 1;
        optopt = 0;
        
        // Create test directories
        fs::create_directories("tests/data/input");
        fs::create_directories("tests/data/output");
        fs::create_directories("tests/data/output/graphs");
        
        // Copy example file to test directory
        fs::copy("../examples/ex1.fasta", "tests/data/input/test.fasta", 
                 fs::copy_options::overwrite_existing);
    }

    void TearDown() override {
        // Reset getopt state after each test
        optind = 1;
        opterr = 1;
        optopt = 0;
        
        // Clean up test files and directories
        std::filesystem::remove_all("tests/data/output");
        std::filesystem::remove("tests/data/input/test.fasta");
    }
};

TEST_F(CliHandlerTest, ShowsHelpMessage) {
    const char* argv[] = {
        "emerald",
        "--help",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;
    
    EXPECT_EXIT(
        CliHandler::parseCommandLine(argc, const_cast<char**>(argv)),
        ::testing::ExitedWithCode(0),
        ".*EMERALD.*Usage.*Required options.*Example.*"
    );
}

TEST_F(CliHandlerTest, ParsesValidCommandLine) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-a", "0.8",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    Config config = CliHandler::parseCommandLine(argc, const_cast<char**>(argv));
    
    EXPECT_EQ(config.input_file, "tests/data/input/test.fasta");
    EXPECT_EQ(config.output_file, "tests/data/output/result.txt");
    EXPECT_FLOAT_EQ(config.alpha, 0.8f);
}

TEST_F(CliHandlerTest, ThrowsOnMissingRequiredOptions) {
    const char* argv[] = {
        "emerald",
        "-a", "0.8",  // Missing required -f and -o options
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    EXPECT_THROW(
        CliHandler::parseCommandLine(argc, const_cast<char**>(argv)),
        std::runtime_error
    );
}

TEST_F(CliHandlerTest, ValidatesAlphaRange) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",  // Use existing test file
        "-o", "tests/data/output/result.txt", // Use proper output path
        "-a", "2.0",  // Invalid alpha value
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    EXPECT_THROW(
        CliHandler::parseCommandLine(argc, const_cast<char**>(argv)),
        std::runtime_error
    );
}

TEST_F(CliHandlerTest, ShowsHelpOnMissingArgs) {
    const char* argv[] = {"emerald", nullptr};
    int argc = 1;

    EXPECT_EXIT({
        try {
            CliHandler::parseCommandLine(argc, const_cast<char**>(argv));
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            exit(1);
        }
    }, ::testing::ExitedWithCode(1),
       ".*EMERALD.*Usage.*Required options.*Input and output files are required.*");
}

TEST_F(CliHandlerTest, ParsesOptionalOutputOptions) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-j", "tests/data/output/report.json",  // JSON output
        "-r", "sequence1",                      // Reference sequence
        "-w", "tests/data/output/graphs",       // Graph output
        "-v",                                   // Verbose mode
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    Config config = CliHandler::parseCommandLine(argc, const_cast<char**>(argv));
    
    EXPECT_EQ(config.json_file, "tests/data/output/report.json");
    EXPECT_EQ(config.reference_protein, "sequence1");
    EXPECT_TRUE(config.draw_graph);
    EXPECT_EQ(config.draw_graph_dir, "tests/data/output/graphs");
    EXPECT_TRUE(config.verbose);
}

TEST_F(CliHandlerTest, ParsesAlignmentParameters) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-d", "5",           // Delta neighborhood
        "-g", "-2",         // Gap cost
        "-e", "-10",        // Start gap
        "-s", "-1",         // Special score
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    Config config = CliHandler::parseCommandLine(argc, const_cast<char**>(argv));
    
    EXPECT_EQ(config.delta, 5);
    EXPECT_EQ(config.gap_cost, -2);
    EXPECT_EQ(config.start_gap, -10);
    EXPECT_EQ(config.special_score, -1);
}

TEST_F(CliHandlerTest, ThrowsOnEmptyInputValue) {
    const char* argv[] = {
        "emerald",
        "-f", "",                               // Empty input path
        "-o", "tests/data/output/result.txt",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    EXPECT_THROW(
        CliHandler::parseCommandLine(argc, const_cast<char**>(argv)),
        std::runtime_error
    );
}

TEST_F(CliHandlerTest, ThrowsOnEmptyOutputValue) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "",                               // Empty output path
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    EXPECT_THROW(
        CliHandler::parseCommandLine(argc, const_cast<char**>(argv)),
        std::runtime_error
    );
}

TEST_F(CliHandlerTest, ThrowsWhenCostMatrixFileMissing) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-c", "tests/data/input/does_not_exist.mat",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    EXPECT_THROW(
        CliHandler::parseCommandLine(argc, const_cast<char**>(argv)),
        std::runtime_error
    );
}

TEST_F(CliHandlerTest, CreatesGraphAndJsonDirectories) {
    namespace fs = std::filesystem;

    // Use non-existent directories; parser should create them
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-w", "tests/data/output/graphs_new",
        "-j", "tests/data/output/reports/report.json",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    // Ensure they don't exist before parsing
    fs::remove_all("tests/data/output/graphs_new");
    fs::remove_all("tests/data/output/reports");

    Config cfg = CliHandler::parseCommandLine(argc, const_cast<char**>(argv));

    EXPECT_TRUE(fs::exists("tests/data/output/graphs_new"));
    EXPECT_TRUE(fs::exists("tests/data/output/reports"));
    EXPECT_EQ(cfg.draw_graph_dir, "tests/data/output/graphs_new");
    EXPECT_EQ(cfg.json_file, "tests/data/output/reports/report.json");
}

TEST_F(CliHandlerTest, SpecialInfSetsIgnoreSpecial) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-s", "INF",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    Config cfg = CliHandler::parseCommandLine(argc, const_cast<char**>(argv));
    EXPECT_TRUE(cfg.ignore_special);
}

TEST_F(CliHandlerTest, WindowMergeFlagSetsTrue) {
    const char* argv[] = {
        "emerald",
        "-f", "tests/data/input/test.fasta",
        "-o", "tests/data/output/result.txt",
        "-m",
        nullptr
    };
    int argc = sizeof(argv) / sizeof(argv[0]) - 1;

    Config cfg = CliHandler::parseCommandLine(argc, const_cast<char**>(argv));
    EXPECT_TRUE(cfg.window_merge);
}