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
#include <cassert>

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi_noop.hpp"
#endif /* USE_MPI */

// MPI params
static int pid, nprocs;

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

    // -- start of spatial partitioning w Morton ordering --

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
    std::vector<std::pair<u64, Star>> keyed;
    keyed.reserve(stars.size());

    for (const auto& s : stars) {
        u64 key = morton2D(s.x, s.y, min_x, min_y, max_x, max_y);
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

    int my_start, my_end;

    if (stars[0].cost == 0) { // first iter
        // naive spatial partition (even split)
        my_start = stars.size() * pid / nprocs;
        my_end   = stars.size() * (pid + 1) / nprocs;
    } else {
        std::vector<int> prefix(stars.size());
        prefix[0] = stars[0].cost;
        for (size_t i = 1; i < stars.size(); i++) {
            prefix[i] = prefix[i-1] + stars[i].cost;
        }

        int total = prefix.back();

        // compute boundaries up front for consistency 
        std::vector<int> boundaries(nprocs + 1);
        boundaries[0] = 0;
        boundaries[nprocs] = stars.size();

        for (int r = 1; r < nprocs; r++) {
            int target = (int)((float)r / nprocs * total);
            // Find first index where prefix >= target
            boundaries[r] = (int)(std::lower_bound(prefix.begin(), prefix.end(), target) 
                                - prefix.begin());
        }

        my_start = boundaries[pid];
        my_end   = boundaries[pid + 1];
    }

    // -- end of stars assignment --

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

    // allgatherv support structures (updated to account for variable counts)
    auto comm_start = chrono::now();
    std::vector<int> counts(nprocs), displs(nprocs);

    int my_count = (my_end - my_start) * sizeof(Star); 
    MPI_Allgather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // for (int vpid = 0; vpid < nprocs; vpid++) {
    //     int node_start = NUM_STARS * vpid / nprocs;
    //     int node_end = NUM_STARS * (vpid + 1) / nprocs;
    //     counts[vpid] = (node_end - node_start) * sizeof(Star);
    // }

    displs[0] = 0;
    for (int vpid = 1; vpid < nprocs; vpid++) {
        displs[vpid] = displs[vpid - 1] + counts[vpid - 1];
    }

    // verify total size matches
    int total_bytes = displs[nprocs-1] + counts[nprocs-1];
    assert(total_bytes == (int)(stars.size() * sizeof(Star)));

    // gather into a separate buffer to avoid aliasing
    std::vector<Star> recv_buf(stars.size());
    MPI_Allgatherv(
        &stars[my_start], my_count, MPI_BYTE,
        recv_buf.data(), counts.data(), displs.data(), MPI_BYTE,
        MPI_COMM_WORLD
    );

    stars = std::move(recv_buf);

    auto comm_end = chrono::now();
    millis comm_time = comm_end - comm_start;

    if (pid == 0) {
        fprintf(stdout, "Qtree took %.01fms, force took %.01fms, comm took %.01fms\n", 
                qtree_time.count(), force_time.count(), comm_time.count());
    }

    destroy_tree(root);
    return total_ct / (my_end - my_start);
}

/**
 * @brief Entry point for single-node MPI implementation
 */
void mpi_single_main(int argc, char* argv[]) {

    fprintf(stderr, "Hello MPI single\n");

#ifndef USE_MPI

    fprintf(stderr, "MPI not enabled!\n");
    exit(-1);

#endif /* USE_MPI */

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

    auto run_start = chrono::now();
    for (int i = 0; i < NUM_ITERS; i++) {
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
    auto run_end = chrono::now();
    millis run_time = run_end - run_start;
    if (pid == 0) display_cleanup();
    MPI_Finalize();

    if (pid == 0) {
        fprintf(stdout, "%d iterations took %.01fms for %.01fms each\n", 
            NUM_ITERS, run_time.count(), run_time.count()/NUM_ITERS);
    }

}