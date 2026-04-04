#include "vector_mpi.hpp"
#include <mpi.h>

/**
 * @brief Send a vector of stars to a specific receiver
 * 
 * @param receiver Receiving node's ID
 * @param stars Stars vector to send
 */
void mpi_send_stars(int receiver, std::vector<Star> stars) {
    int length = stars.size();
    MPI_Send(&length, 1, MPI_INT, receiver, 0, MPI_COMM_WORLD);

    if (length > 0) {
        MPI_Send(stars.data(), sizeof(Star)*stars.size(), MPI_BYTE, receiver, 1, MPI_COMM_WORLD);
    }
}

/**
 * @brief Receive a stars vector
 * 
 * @param sender Sender node's ID
 * @return Received vector
 */
std::vector<Star> mpi_recv_stars(int sender) {
    int length;
    MPI_Recv(&length, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    std::vector<Star> my_stars(length);
    if (length > 0) {
        MPI_Recv(my_stars.data(), sizeof(Star)*my_stars.size(), MPI_BYTE, sender, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    return my_stars;
}
