#pragma once
#include "../core/config.h"
#include <string>

class CliHandler {
public:
    static Config parseCommandLine(int argc, char** argv);
    
private:
    static void printUsage(const char* programName);
    static void validateArguments(const Config& config);
};