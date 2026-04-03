#include "compact_defines.h"
#include "serial.hpp"
#include "display.hpp"

/**
 * @brief Entry point for Heterogeneous Distributed Barnes-Hut
 */
int main(int argc, char* argv[]) {

    if (argc != 2) {
        fprintf(stderr, "Need argument\n");
        return -1;
    }

    display_init();

    if (strcmp(argv[1], "serial") == 0) {
        serial_main();
        return 0;
    } else if (strcmp(argv[1], "ghc") == 0) {
        fprintf(stderr, "Running on GHC\n");
        return 0;
    } else if (strcmp(argv[1], "ece") == 0) {
        fprintf(stderr, "Running on ECE\n");
        return 0;
    } else if (strcmp(argv[1], "linux") == 0) {
        fprintf(stderr, "Running on linux.andrew\n");
        return 0;
    } else if (strcmp(argv[1], "hh") == 0) {
        fprintf(stderr, "Running on HH\n");
        return 0;
    }

    fprintf(stderr, "Need argument\n");
    return -1;
}
