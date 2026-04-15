/**
 * @file minimpi.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Header for Mini-MPI
 */

#pragma once

#include "compact_defines.h"
#include <vector>

int mmpi_getpid(void);
int mmpi_getnodes(void);

void mmpi_init(int init_pid, int init_node_ct);
void mmpi_finalize(void);
void mmpi_send(int receiver, void* buf, int length);
void mmpi_recv(int sender, void* buf, int length);
void mmpi_bcast(int sender, void* buf, int length);
void mmpi_barrier(void);
void mmpi_sync(void* buf, int length, int contrib_bytes);
void mmpi_syncv(void* buf, int length, int contribs[], int displs[]);

/**
 * @brief Send a vector to a specific receiver
 * 
 * @param receiver Receiving node's ID
 * @param vector Vector to send
 */
template <typename T>
void mmpi_send_vec(int receiver, const std::vector<T>& vec) {
    int length = vec.size();

    // send length first
    mmpi_send(receiver, &length, sizeof(int));

    // send data
    if (length > 0) {
        mmpi_send(receiver, (void*)vec.data(), sizeof(T) * length);
    }
}

/**
 * @brief Receive a vector
 * 
 * @param sender Sender node's ID
 * @return Received vector
 */
template <typename T>
std::vector<T> mmpi_recv_vec(int sender) {
    int length;

    // receive length
    mmpi_recv(sender, &length, sizeof(int));

    std::vector<T> vec(length);

    // receive data
    if (length > 0) {
        mmpi_recv(sender, vec.data(), sizeof(T) * length);
    }

    return vec;
}
