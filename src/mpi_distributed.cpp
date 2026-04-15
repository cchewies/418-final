/**
 * @file mpi_distributed.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * @author Vanessa Lam (yatheil@andrew.cmu.edu)
 * 
 * Multi-node Barnes-Hut
 */

#include "mpi_distributed.hpp"
#include "compact_defines.h"
#include "minimpi.hpp"
#include "simulation_config.hpp"
#include "display.hpp"
#include "quadtree.hpp"
#include <cmath>
#include <algorithm>

// MPI params
static int pid, nprocs;

/**
 * @brief Quadtree simulation step
 * 
 * @param stars Vector of stars
 * @return average # of stars visited
 */
static int mpi_iterate_simulation(std::vector<Star> &stars) {

    // -- start of spatial partitioning w Morton ordering --

    auto assign_start = chrono::now();

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

    // Unbalanced
    // get node's allocation (pure naive)
    my_start = NUM_STARS * (pid) / nprocs;
    my_end = NUM_STARS * (pid+1) / nprocs;

    // Load balanced
    // if (stars[0].cost == 0) { // first iter
    //     // naive allocation
    //     my_start = stars.size() * (pid) / nprocs;
    //     my_end = stars.size() * (pid+1) / nprocs;
    // } else {
    //     // prefix-sum load balancing

    //     std::vector<int> prefix(stars.size());
    //     prefix[0] = stars[0].cost;
    //     for (size_t i = 1; i < stars.size(); i++) {
    //         prefix[i] = prefix[i-1] + stars[i].cost;
    //     }

    //     int total = prefix.back();

    //     // compute boundaries up front for consistency 
    //     std::vector<int> boundaries(nprocs + 1);
    //     boundaries[0] = 0;
    //     boundaries[nprocs] = stars.size();

    //     for (int r = 1; r < nprocs; r++) {
    //         int target = (int)((float)r / nprocs * total);
    //         // Find first index where prefix >= target
    //         boundaries[r] = (int)(std::lower_bound(prefix.begin(), prefix.end(), target) 
    //                             - prefix.begin());
    //     }

    //     my_start = boundaries[pid];
    //     my_end   = boundaries[pid + 1];
    // }

    // -- end of stars assignment --
    auto assign_end = chrono::now();
    millis assign_time = assign_end - assign_start;

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

    // Load balanced
    // int my_count = (my_end - my_start) * sizeof(Star);
    // counts[pid] = my_count;
    // mmpi_sync(counts.data(), counts.size() * sizeof(int), sizeof(int));

    // Unbalanced
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
    mmpi_syncv(stars.data(), stars.size() * sizeof(Star), counts.data(), displs.data());

    auto comm_end = chrono::now();
    millis comm_time = comm_end - comm_start;

    fprintf(stdout, "Assign took %.01fms, qtree took %.01fms, force took %.01fms, comm took %.01fms\n", 
            assign_time.count(), qtree_time.count(), force_time.count(), comm_time.count());

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
    pid = mmpi_getpid();
    nprocs = mmpi_getnodes();

    if (pid == 0) display_init();

    std::vector<Star> stars;
    if (pid == 0) {
        stars = generate_galaxy();
        for (int receiver = 1; receiver < nprocs; receiver++) {
            mmpi_send_vec<Star>(receiver, stars);
        }
    } else {
        stars = mmpi_recv_vec<Star>(0);
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

        fprintf(stdout, "Iteration took %.01fms, avg stars: %d/%ld\n", 
            frame_time.count(), avg_stars, stars.size());

        bool quit = false;
        if (pid == 0) {
            quit = check_quit();
        }
        mmpi_bcast(0, &quit, sizeof(bool)); // TODO: this is probably also adding some latency
        if (quit) break;
    }
    auto run_end = chrono::now();
    millis run_time = run_end - run_start;
    if (pid == 0) display_cleanup();

    mmpi_finalize();

    if (pid == 0) {
        fprintf(stdout, "%d iterations took %.01fms for %.01fms each\n", 
            NUM_ITERS, run_time.count(), run_time.count()/NUM_ITERS);
    }
}
