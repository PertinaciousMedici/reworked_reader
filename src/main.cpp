#include "parser.h"
#include <iostream>

int main(int argc, char** argv) {
    if (argc <= 1) {
        std::cout << "Usage: " << argv[0] << " <filepath>" << std::endl;
        return 1;
    }

    try {
        auto instance = parseTSP(argv[1]);
        size_t dimension = instance.dimension;

        std::cout << dimension << std::endl;
        std::cout << "Distance Matrix: " << std::endl;
        printMatrix(instance);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
        return 1;
    }

    return 0;
}