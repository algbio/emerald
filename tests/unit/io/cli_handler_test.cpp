#include <gtest/gtest.h>
#include "../../../src/io/cli_handler.h"
#include "../../../src/core/config.h"
#include <string>
#include <cstdio>
#include <filesystem>

class CliHandlerTest : public ::testing::Test {
protected:
    void SetUp() override {
        namespace fs = std::filesystem;
        
        // Create test directories
        fs::create_directories("tests/data/input");
        fs::create_directories("tests/data/output");
        
        // Copy example file to test directory
        fs::copy("../examples/ex1.fasta", "tests/data/input/test.fasta", 
                 fs::copy_options::overwrite_existing);
    }

    void TearDown() override {
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
        "-f", "input.fasta",
        "-o", "output.txt",
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