#include "mpi_distributed.hpp"
#include "compact_defines.h"
#include "minimpi.hpp"
#include "simulation_config.hpp"
#include "display.hpp"
#include "quadtree.hpp"
#include <cmath>

/**
 * @brief Barnes-hut force calculation
 * 
 * @param s Star to calculate total force for
 * @param node Node to add forces to
 * @param fx Reference to x force accumulator
 * @param fy Reference to y force accumulator
 * @return # of stars visited
 */
static int compute_force(Star& s, QNode* node, float& fx, float& fy) {

    int star_count = 0;

    // Node doesnt exist, no stars, no children, or star is current one
    if (!node || node->mass == 0 || (node->s == &s && !node->nw)) {
        return star_count;
    }

    float dx = node->com_x - s.x;
    float dy = node->com_y - s.y;
    float dist2 = dx*dx + dy*dy + EPS2;
    float dist = std::sqrt(dist2);

    if (!node->nw || node->side_len / dist < THETA) {
        // no children
        float force = G * s.mass * node->mass / dist2;
        fx += force * dx / dist;
        fy += force * dy / dist;
        star_count++;
    } else {
        // children exist
        star_count += compute_force(s, node->nw, fx, fy);
        star_count += compute_force(s, node->ne, fx, fy);
        star_count += compute_force(s, node->sw, fx, fy);
        star_count += compute_force(s, node->se, fx, fy);
    }

    return star_count;
}

/**
 * @brief Quadtree simulation step
 * 
 * @param stars Vector of stars
 * @return average # of stars visited
 */
static int mpi_iterate_simulation(std::vector<Star> &stars) {

    auto qtree_start = chrono::now();
    QNode* root = build_qtree(stars);
    auto qtree_end = chrono::now();
    millis qtree_time = qtree_end - qtree_start;

    int total_ct = 0;

    // get node's allocation
    int my_start = NUM_STARS * (mmpi_getpid()) / mmpi_getnodes();
    int my_end = NUM_STARS * (mmpi_getpid()+1) / mmpi_getnodes();

    // update velocities
    auto force_start = chrono::now();
    for (int i = my_start; i < my_end; i++) {
        Star& s = stars[i];
        float fx = 0, fy = 0;
        int star_ct = compute_force(s, root, fx, fy);
        total_ct += star_ct;

        s.vx += fx / s.mass * DT;
        s.vy += fy / s.mass * DT;
    }
    auto force_end = chrono::now();
    millis force_time = force_end - force_start;

    // update positions
    for (int i = my_start; i < my_end; i++) {
        Star& s = stars[i];
        s.x += s.vx * DT;
        s.y += s.vy * DT;
    }

    // allgatherv support structures
    std::vector<int> counts(mmpi_getnodes()), displs(mmpi_getnodes());
    for (int vpid = 0; vpid < mmpi_getnodes(); vpid++) {
        int node_start = NUM_STARS * vpid / mmpi_getnodes();
        int node_end = NUM_STARS * (vpid + 1) / mmpi_getnodes();
        counts[vpid] = (node_end - node_start) * sizeof(Star);
    }

    displs[0] = 0;
    for (int vpid = 1; vpid < mmpi_getnodes(); vpid++) {
        displs[vpid] = displs[vpid - 1] + counts[vpid - 1];
    }

    // gather updated stars to all ranks
    mmpi_syncv(stars.data(), stars.size() * sizeof(Star), counts.data(), displs.data());

    fprintf(stdout, "Qtree took %.01fms, force took %.01fms\n", 
            qtree_time.count(), force_time.count());

    destroy_tree(root);
    return total_ct / (my_end - my_start);
}

/**
 * @brief Entry point for distributed mpi implementation
 */
void mpi_distributed_main(int argc, char* argv[]) {

    printf("Hello MPI distributed\n");

    if (argc != 4) {
        fprintf(stderr, "Need 3 arguments for distributed\n");
        exit(-1);
    }

    int init_pid = atoi(argv[2]);
    int init_node_ct = atoi(argv[3]);

    mmpi_init(init_pid, init_node_ct);

    if (mmpi_getpid() == 0) display_init();

    std::vector<Star> stars;
    if (mmpi_getpid() == 0) {
        stars = generate_galaxy();
        for (int receiver = 1; receiver < mmpi_getnodes(); receiver++) {
            mmpi_send_vec<Star>(receiver, stars);
        }
    } else {
        stars = mmpi_recv_vec<Star>(0);
    }

    while (1) {
        if (mmpi_getpid() == 0) {
            display_render(stars);
        }

        auto start = chrono::now();

        int avg_stars = mpi_iterate_simulation(stars);
        
        auto end = chrono::now();
        millis frame_time = end - start;

        fprintf(stdout, "Iteration took %.01fms, avg stars: %d/%ld\n", 
            frame_time.count(), avg_stars, stars.size());

        bool quit = false;
        if (mmpi_getpid() == 0) {
            quit = check_quit();
        }
        mmpi_bcast(0, &quit, sizeof(bool));
        if (quit) break;
    }
    if (mmpi_getpid() == 0) display_cleanup();

    mmpi_finalize();
}
