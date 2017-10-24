#include <fstream>

#include "call.hpp"

int main(int argc, char** argv) {
    if (argc > 1) {
        auto in = std::ifstream(argv[1]);
        callBayes(in);
    }
}