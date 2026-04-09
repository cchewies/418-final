#include "minimpi.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <limits.h>
#include <cassert>

#define PORT 5000
#define MAX_NODES 16
#define COORDINATOR_HOSTNAME ("ghc82.ghc.andrew.cmu.edu")

static char hostname_buf[MAX_NODES][HOST_NAME_MAX];
static int socks[MAX_NODES] = {-1};

static int pid;
static int num_nodes;

int mmpi_getpid(void) {
    return pid;
}

int mmpi_getnodes(void) {
    return num_nodes;
}

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

    for (int i = 1; i < num_nodes; i++) {

        listen(coord_fd, MAX_NODES);
        printf("Listening on port %d\n", PORT);

        // Accept connections from workers
        int sock = accept(coord_fd, (struct sockaddr*)&addr, &addrlen);
        if (sock < 0) {
            perror("accept"); 
            exit(1);
        }
        
        // Receive message
        int recv_pid;
        recv(sock, &recv_pid, sizeof(recv_pid), 0);

        assert(recv_pid < num_nodes);
        socks[recv_pid] = sock;

        printf("Worker %d (%d/%d) connected!\n", recv_pid, i, num_nodes-1);

        recv(sock, hostname_buf[recv_pid], sizeof(hostname_buf[recv_pid]), 0);
        printf("Received from worker: %s\n", hostname_buf[recv_pid]);
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
    // Get coordinator ip
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
    send(socks[0], &pid, sizeof(pid), 0);
    send(socks[0], hostname_buf[pid], strlen(hostname_buf[pid])+1, 0);

    // Receive message
    recv(socks[0], hostname_buf, sizeof(hostname_buf), 0);
    printf("Received from coordinator: %s\n", hostname_buf[0]);
    printf("Received from coordinator: %s\n", hostname_buf[1]);
    printf("Received from coordinator: %s\n", hostname_buf[2]);
    printf("Received from coordinator: %s\n", hostname_buf[3]);

    // Workers can now open pairwise connections
    int listen_fd = socket(AF_INET, SOCK_STREAM, 0);
    if (listen_fd < 0) { perror("socket"); exit(1); }

    // Open a listening connection to all pids higher
    for (int serve = pid+1; serve < num_nodes; serve++) {
        printf("%d listening for %d: %s\n", pid, serve, hostname_buf[serve]);

        struct sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_addr.s_addr = INADDR_ANY;
        addr.sin_port = htons(PORT + pid);
        bind(listen_fd, (struct sockaddr*)&addr, sizeof(addr));
        listen(listen_fd, 1);
        if (listen_fd < 1) { perror("listen"); exit(1); }

        // store connection to worker
        int conn_fd = accept(listen_fd, NULL, NULL);
        if (conn_fd < 1) { perror("accept"); exit(1); }
        socks[serve] = conn_fd;
    }
    
    // the listening socket can be closed now
    close(listen_fd);

    // Ping all pids lower
    for (int ping = pid-1; ping >= 1; ping--) {
        printf("%d pinging %d: %s\n", pid, ping, hostname_buf[ping]);
        int sock = socket(AF_INET, SOCK_STREAM, 0);
        if (sock < 1) { perror("sock"); exit(1); }

        struct sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_port = htons(PORT + ping);

        // resolve hostname of worker (from coordinator info)
        struct hostent *he = gethostbyname(hostname_buf[ping]);
        memcpy(&addr.sin_addr, he->h_addr_list[0], he->h_length);

        // retry until connected
        while (connect(sock, (struct sockaddr*)&addr, sizeof(addr)) < 0) {
            perror("connect");
            usleep(100000);
        }
        printf("Ping success\n");

        socks[ping] = sock;
    }
    
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
 * @brief Finalize mmpi resources
 */
void mmpi_finalize(void) {
    for (int i = 0; i < num_nodes; i++) {
        if (i == pid) {
            continue;
        }
        
        assert(socks[i] >= 0);

        close(socks[i]);
        socks[i] = -1;
    }
    printf("%d: all sockets closed\n", pid);
}

/**
 * @brief Send a message to a specified node
 * 
 * @param receiver Receiving node's pid
 * @param buf Buffer of data to send
 * @param length Length of buf in bytes
 */
void mmpi_send(int receiver, void* buf, int length) {
    assert(receiver != pid);
    assert(receiver < num_nodes);

    int sock = socks[receiver];
    int sent = 0;

    while (sent < length) {
        int n = send(sock, (char*)buf + sent, length - sent, 0);
        if (n <= 0) {
            perror("send");
            exit(1);
        }
        sent += n;
    }
}

/**
 * @brief Receive a message from a specified node
 * 
 * @param sender Sender node's pid
 * @param buf Buffer to receive into
 * @param length Length of message to receive in bytes
 */
void mmpi_recv(int sender, void* buf, int length) {
    assert(sender != pid);
    assert(sender < num_nodes);

    int sock = socks[sender];
    int received = 0;

    while (received < length) {
        int n = recv(sock, (char*)buf + received, length - received, 0);
        if (n <= 0) {
            perror("recv");
            exit(1);
        }
        received += n;
    }
}

/**
 * @brief Broadcast or receive a message (all nodes)
 * 
 * @param sender Sender of broadcast
 * @param buf Buffer to send/receive
 * @param length Length of buffer in bytes
 */
void mmpi_bcast(int sender, void* buf, int length) {

    if (pid == sender) {
        // Sender needs to send data to all nodes

        for (int receiver = 0; receiver < num_nodes; receiver++) {
            if (pid == receiver) {
                continue;
            }

            mmpi_send(receiver, buf, length);
        }

    } else {
        // Receivers wait for data
        mmpi_recv(sender, buf, length);
    }

}

/**
 * @brief Gathers data from all nodes and delivers it to all
 *        Each process may contribute a different amount of data
 * 
 * @param buf Buffer to send from / recv into
 * @param length Length of buffer in bytes
 * @param contribs Length of each pid's contributions 
 * @param displs Displacement of each pid's region
 */
void mmpi_syncv(void* buf, int length, int contribs[], int displs[]) {
    // TODO: temporary naive implementation
    if (pid == 0) {
        for (int i = 0; i < num_nodes; i++) {
            if (i == pid) {
                continue;
            }

            mmpi_recv(i, (char*)buf + displs[i], contribs[i]);
        }
    } else {
        mmpi_send(0, (char*)buf + displs[pid], contribs[pid]);
    }

    mmpi_bcast(0, buf, length);
}
