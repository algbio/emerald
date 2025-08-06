#include <iostream>
#include "core/application.h"
#include "io/cli_handler.h"

int main(int argc, char** argv) {
    try {
        Config config = CliHandler::parseCommandLine(argc, argv);
        Application app(config);
        return app.run();
    } 
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
