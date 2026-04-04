/**
 * @file entry.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Entry point for all of the various Barnes-Hut implementations.
 * 
 * naive: O(N^2), single-thread, single-node
 * serial: O(NlogN), single-thread, single-node
 * mpi: O(NlogN), single-node MPI
 * distributed: O(NlogN), distributed MPI
 */

#include "compact_defines.h"
#include "serial.hpp"
#include "mpi_single.hpp"
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

    if (strcmp(argv[1], "naive") == 0) {
        naive_main(argc, argv);
        return 0;
    } else if (strcmp(argv[1], "serial") == 0) {
        serial_main(argc, argv);
        return 0;
    } else if (strcmp(argv[1], "mpi") == 0) {
        mpi_main(argc, argv);
        return 0;
    } else if (strcmp(argv[1], "distributed") == 0) {
        fprintf(stderr, "distributed implementation\n");
        fprintf(stderr, "not implemented\n");
        return -1;
    }

    fprintf(stderr, "Need argument\n");
    return -1;
}
