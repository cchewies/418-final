#include "minimpi.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <limits.h>

#define PORT 5000
#define MAX_NODES 4
#define COORDINATOR_HOSTNAME ("ghc82.ghc.andrew.cmu.edu")

static char hostname_buf[MAX_NODES][HOST_NAME_MAX];
static int socks[MAX_NODES] = {-1};

static int pid;
static int num_nodes;

/**
 * @brief PID 0 is the coordinator
 */
static bool is_coordinator(void) {
    return pid == 0;
}

/**
 * @brief Coordinator should send information for pairwise connections
 */
static void mmpi_init_coordinator(void) {

    memcpy(hostname_buf[0], COORDINATOR_HOSTNAME, sizeof(COORDINATOR_HOSTNAME));

    // Create TCP socket
    int coord_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (coord_fd < 0) { 
        perror("socket");
        exit(1); 
    }

    // Bind to port
    struct sockaddr_in addr;
    socklen_t addrlen = sizeof(addr);
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = INADDR_ANY;
    addr.sin_port = htons(PORT);
    if (bind(coord_fd, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
        perror("bind");
        exit(1);
    }

    for (int worker = 1; worker < num_nodes; worker++) {

        listen(coord_fd, MAX_NODES);
        fprintf(stderr, "Listening on port %d\n", PORT);

        // Accept connections from workers
        socks[worker] = accept(coord_fd, (struct sockaddr*)&addr, &addrlen);
        if (socks[worker] < 0) {
            perror("accept"); 
            exit(1);
        }

        fprintf(stderr, "Worker %d/%d connected!\n", worker, num_nodes-1);

        // Receive message
        recv(socks[worker], hostname_buf[worker], sizeof(hostname_buf[worker]), 0);
        printf("Received from worker: %s\n", hostname_buf[worker]);
    }

    for (int worker = 1; worker < num_nodes; worker++) {
        // Send messages
        send(socks[worker], hostname_buf, sizeof(hostname_buf), 0);
    }
    
    close(coord_fd);
}

/**
 * @brief Workers should send info to coordinator, then receive connections
 */
static void mmpi_init_worker(void) {
    // Resolve hostname to IP
    struct hostent* he = gethostbyname(COORDINATOR_HOSTNAME);
    if (he == NULL) {
        perror("gethostbyname");
        exit(1);
    }

    socks[0] = socket(AF_INET, SOCK_STREAM, 0);
    if (socks[0] < 0) {
        perror("socket"); 
        exit(1);
    }

    struct sockaddr_in server_addr;
    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(PORT);
    memcpy(&server_addr.sin_addr, he->h_addr_list[0], he->h_length);

    // Connect to coordinator
    if (connect(socks[0], (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        perror("connect");
        exit(1);
    }
    
    // Send message
    gethostname(hostname_buf[pid], sizeof(hostname_buf[pid]));
    send(socks[0], hostname_buf[pid], strlen(hostname_buf[pid])+1, 0);

    // Receive message
    recv(socks[0], hostname_buf, sizeof(hostname_buf), 0);
    printf("Received from coordinator: %s\n", hostname_buf[0]);
    printf("Received from coordinator: %s\n", hostname_buf[1]);
    printf("Received from coordinator: %s\n", hostname_buf[2]);
    printf("Received from coordinator: %s\n", hostname_buf[3]);
}
    

/**
 * @brief Initialize minimpi
 * 
 * @param init_pid PID supplied from cli args
 * @param init_node_ct Number of nodes in this instance
 */
void mmpi_init(int init_pid, int init_node_ct) {

    pid = init_pid;
    num_nodes = init_node_ct;

    printf("my pid %d, with %d nodes\n", pid, num_nodes);

    if (is_coordinator()) {
        mmpi_init_coordinator();
    } else {
        mmpi_init_worker();
    }

}

/**
 * @brief Send a message to a specified node
 * 
 * @param receiver Receiving node's pid
 * @param length Length of buf in bytes
 * @param buf Buffer of data to send
 */
void mmpi_send(int receiver, int length, void* buf) {

}

/**
 * @brief Receive a message from a specified node
 * 
 * @param sender Sender node's pid
 * @param length Length of message to receive
 * @param buf Buffer to receive into
 */
void mmpi_recv(int sender, int length, void* buf) {

}
