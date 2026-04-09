/**
 * @file mpi_single.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * @author Vanessa Lam (yatheil@andrew.cmu.edu)
 * 
 * Single-node MPI Barnes-Hut
 * 
 * mpi: O(NlogN), single-node MPI
 */

#include "compact_defines.h"
#include "simulation_config.hpp"
#include "mpi_single.hpp"
#include "vector_mpi.hpp"
#include "display.hpp"
#include "quadtree.hpp"
#include <cmath>
#include <cstdint>
#include <algorithm>

#ifdef USE_MPI
#include <mpi.h>

// MPI params
static int pid, nprocs;

// expand 32-bit int into 64-bit w interleaved zeros
static inline uint64_t expand_bits(uint32_t x) {
    uint64_t v = x;
    v = (v | (v << 16)) & 0x0000FFFF0000FFFFULL;
    v = (v | (v << 8))  & 0x00FF00FF00FF00FFULL;
    v = (v | (v << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    v = (v | (v << 2))  & 0x3333333333333333ULL;
    v = (v | (v << 1))  & 0x5555555555555555ULL;
    return v;
}

// 2D Morton encoding
static inline uint64_t morton2D(float x, float y,
                                float min_x, float min_y,
                                float max_x, float max_y) {
    float nx = (x - min_x) / (max_x - min_x); // normalize to [0, 1]
    float ny = (y - min_y) / (max_y - min_y);
    nx = std::min(1.0f, std::max(0.0f, nx));
    ny = std::min(1.0f, std::max(0.0f, ny));
    uint32_t ix = (uint32_t)(nx * ((1 << 21) - 1));
    uint32_t iy = (uint32_t)(ny * ((1 << 21) - 1));
    return (expand_bits(ix) << 1) | expand_bits(iy);
}

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

    // naive allocation
    // int my_start = NUM_STARS * (pid) / nprocs;
    // int my_end = NUM_STARS * (pid+1) / nprocs;

    // start of spatial partitioning w Morton ordering 

    // compute bounding box
    float min_x = stars[0].x, max_x = stars[0].x;
    float min_y = stars[0].y, max_y = stars[0].y;

    for (const auto& s : stars) {
        min_x = std::min(min_x, s.x);
        max_x = std::max(max_x, s.x);
        min_y = std::min(min_y, s.y);
        max_y = std::max(max_y, s.y);
    }

    // compute keys
    std::vector<std::pair<uint64_t, Star>> keyed;
    keyed.reserve(stars.size());

    for (const auto& s : stars) {
        uint64_t key = morton2D(s.x, s.y, min_x, min_y, max_x, max_y);
        keyed.emplace_back(key, s);
    }

    // sort by Morton key
    std::sort(keyed.begin(), keyed.end(),
        [](const auto& a, const auto& b) {
            return a.first < b.first;
        }
    );

    // write back sorted stars
    for (size_t i = 0; i < stars.size(); i++) {
        stars[i] = keyed[i].second;
    }

    // compute spatial partition
    int my_start = stars.size() * pid / nprocs;
    int my_end   = stars.size() * (pid + 1) / nprocs;

    // end of stars assignment 

    auto qtree_start = chrono::now();
    QNode* root = build_qtree(stars);
    auto qtree_end = chrono::now();
    millis qtree_time = qtree_end - qtree_start;

    int total_ct = 0;

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
    auto comm_start = chrono::now();
    std::vector<int> counts(nprocs), displs(nprocs);
    for (int vpid = 0; vpid < nprocs; vpid++) {
        int node_start = NUM_STARS * vpid / nprocs;
        int node_end = NUM_STARS * (vpid + 1) / nprocs;
        counts[vpid] = (node_end - node_start) * sizeof(Star);
    }

    displs[0] = 0;
    for (int vpid = 1; vpid < nprocs; vpid++) {
        displs[vpid] = displs[vpid - 1] + counts[vpid - 1];
    }

    // gather updated stars to all ranks
    MPI_Allgatherv(
        &stars[my_start], (my_end - my_start) * sizeof(Star), MPI_BYTE,
        stars.data(), counts.data(), displs.data(), MPI_BYTE,
        MPI_COMM_WORLD
    );
    auto comm_end = chrono::now();
    millis comm_time = comm_end - comm_start;

    if (pid == 0) {
        fprintf(stdout, "Qtree took %.01fms, force took %.01fms, comm took %.01fms\n", 
                qtree_time.count(), force_time.count(), comm_time.count());
    }

    destroy_tree(root);
    return total_ct / (my_end - my_start);
}

#endif /* USE_MPI */

/**
 * @brief Entry point for single-node MPI implementation
 */
void mpi_single_main(int argc, char* argv[]) {

    fprintf(stderr, "Hello MPI single\n");

#ifdef USE_MPI
    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (pid == 0) display_init();

    std::vector<Star> stars;
    if (pid == 0) {
        stars = generate_galaxy();
        for (int receiver = 1; receiver < nprocs; receiver++) {
            mpi_send_stars(receiver, stars);
        }
    } else {
        stars = mpi_recv_stars(0);
    }

    while (1) {
        if (pid == 0) {
            display_render(stars);
        }

        auto start = chrono::now();

        int avg_stars = mpi_iterate_simulation(stars);
        
        auto end = chrono::now();
        millis frame_time = end - start;

        if (pid == 0) {
            fprintf(stdout, "Iteration took %.01fms, avg stars: %d/%ld\n", 
                frame_time.count(), avg_stars, stars.size());
        }

        bool quit = false;
        if (pid == 0) {
            quit = check_quit();
        }
        MPI_Bcast(&quit, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        if (quit) break;
    }
    if (pid == 0) display_cleanup();
    MPI_Finalize();

#else

    fprintf(stderr, "MPI not enabled!\n");
    exit(-1);

#endif /* USE_MPI */
}