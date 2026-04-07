#include "mpi_distributed.hpp"
#include "compact_defines.h"
#include "minimpi.hpp"

void mpi_distributed_main(int argc, char* argv[]) {

    fprintf(stderr, "Hello MPI distributed\n");

    if (argc != 4) {
        fprintf(stderr, "Need 3 arguments for distributed\n");
        exit(-1);
    }

    int init_pid = atoi(argv[2]);
    int init_node_ct = atoi(argv[3]);

    mmpi_init(init_pid, init_node_ct);
}
