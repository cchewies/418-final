/**
 * @file mpi_parallel.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * @author Vanessa Lam (yatheil@andrew.cmu.edu)
 * 
 * Single-node MPI Barnes-Hut
 * Multi-node MMPI Barnes-Hut
 */

#include "mpi_parallel.hpp"
#include "compact_defines.h"
#include "minimpi.hpp"
#include "simulation_config.hpp"
#include "display.hpp"
#include "quadtree.hpp"
#include "vector_mpi.hpp"
#include <cmath>
#include <algorithm>

#ifdef USE_MPI
#include <mpi.h>
#else
#include "mpi_noop.hpp"
#endif /* USE_MPI */

// MPI params
static int pid, nprocs;
static bool use_mmpi;


std::vector<std::pair<u64, Star>> keyed;

/**
 * @brief Quadtree simulation step
 * 
 * @param stars Vector of stars
 * @param positions Preallocated positions vector
 * @return average # of stars visited
 */
static int mpi_iterate_simulation(std::vector<Star> &stars,
                                  std::vector<StarPos> &positions) {

    // -- start of spatial partitioning w Morton ordering --

    auto assign_start = chrono::now();

    // compute bounding box
    float min_x = stars[0].pos.x, max_x = stars[0].pos.x;
    float min_y = stars[0].pos.y, max_y = stars[0].pos.y;

    for (const auto& s : stars) {
        min_x = std::min(min_x, s.pos.x);
        max_x = std::max(max_x, s.pos.x);
        min_y = std::min(min_y, s.pos.y);
        max_y = std::max(max_y, s.pos.y);
    }

    // compute keys
    keyed.reserve(stars.size());

    for (const auto& s : stars) {
        u64 key = morton2D(s.pos.x, s.pos.y, min_x, min_y, max_x, max_y);
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
    // update positions and put into aux struct
    for (int i = my_start; i < my_end; i++) {
        Star& s = stars[i];
        s.pos.x += s.vx * DT;
        s.pos.y += s.vy * DT;
        positions[i+pid+1] = s.pos;
    }
    auto force_end = chrono::now();
    millis force_time = force_end - force_start;
    millis compute_time = force_end - assign_start;

    // allgatherv support structures
    auto comm_start = chrono::now();
    std::vector<int> counts(nprocs), displs(nprocs);
    int my_count = (my_end - my_start + 1) * sizeof(StarPos);       // +1 for timing

    // Unbalanced
    for (int vpid = 0; vpid < nprocs; vpid++) {
        int node_start = NUM_STARS * vpid / nprocs;
        int node_end = NUM_STARS * (vpid + 1) / nprocs;
        counts[vpid] = (node_end - node_start + 1) * sizeof(StarPos);   // +1 for timing
    }

    // Load balanced
    // counts[pid] = my_count;
    // if (use_mmpi) {
    //     mmpi_sync(counts.data(), counts.size() * sizeof(int), sizeof(int));
    // } else {
    //     MPI_Allgather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    // }

    displs[0] = 0;
    for (int vpid = 1; vpid < nprocs; vpid++) {
        displs[vpid] = displs[vpid - 1] + counts[vpid - 1];
    }
    positions[displs[pid]/((int)sizeof(StarPos))].runtime = compute_time.count();

    // gather updated stars to all ranks
    if (use_mmpi) {
        mmpi_syncv(positions.data(), positions.size() * sizeof(StarPos), counts.data(), displs.data());
    } else {
        // gather into a separate buffer to avoid aliasing
        std::vector<StarPos> recv_buf(positions.size());
        MPI_Allgatherv(
            &positions[my_start], my_count, MPI_BYTE,
            recv_buf.data(), counts.data(), displs.data(), MPI_BYTE,
            MPI_COMM_WORLD
        );

        positions = std::move(recv_buf);
    }

    int star_idx = 0;
    int pos_idx = 0;
    for (int vpid = 0; vpid < nprocs; vpid++) {
        pos_idx++;

        int num_elems = counts[vpid] / ((int)sizeof(StarPos)) - 1;
        for (int i = 0; i < num_elems; i++) {
            stars[star_idx].pos = positions[pos_idx];
            // printf("%f\n", stars[star_idx].pos.x);
            star_idx++;
            pos_idx++;
        }
    }
    // for (int i = 0; i < NUM_STARS; i++) {
    //     stars[i].pos = positions[i];
    // }

    printf("I took %.01fms\n", compute_time.count());
    std::vector<double> compute_times(nprocs);

    // for (auto p : positions) {
    //     printf("%f %f\n", p.x, p.runtime);
    // }

    for (int i = 0; i < nprocs; i++) {
        compute_times[i] = positions[displs[i]/((int)sizeof(StarPos))].runtime;
        printf("pid %d took %.01fms\n", i, compute_times[i]);
    }
    // mmpi_sync(compute_times.data(), compute_times.size() * sizeof(double), sizeof(double));

    auto comm_end = chrono::now();
    millis comm_time = comm_end - comm_start;

    if (use_mmpi || pid == 0) {
        fprintf(stdout, "Assign took %.01fms, qtree took %.01fms, force took %.01fms, comm took %.01fms\n", 
                assign_time.count(), qtree_time.count(), force_time.count(), comm_time.count());
    }

    destroy_tree(root);
    return total_ct / (my_end - my_start);
}

/**
 * @brief Run the simulation
 */
static void mpi_run_simulation(void) {
    if (pid == 0) display_init();

    std::vector<Star> stars;
    if (pid == 0) {
        stars = generate_galaxy();
        for (int receiver = 1; receiver < nprocs; receiver++) {
            if (use_mmpi) {
                mmpi_send_vec<Star>(receiver, stars);
            } else {
                mpi_send_stars(receiver, stars);
            }
        }
    } else {
        if (use_mmpi) {
            stars = mmpi_recv_vec<Star>(0);
        } else {
            stars = mpi_recv_stars(0);
        }
    }

    std::vector<StarPos> positions(NUM_STARS+nprocs);
    positions.resize(NUM_STARS+nprocs);

    auto run_start = chrono::now();
    for (int i = 0; i < NUM_ITERS; i++) {
        if (pid == 0) {
            display_render(stars);
        }

        auto start = chrono::now();

        int avg_stars = mpi_iterate_simulation(stars, positions);
        
        auto end = chrono::now();
        millis frame_time = end - start;

        fprintf(stdout, "Iteration took %.01fms, avg stars: %d/%ld\n", 
            frame_time.count(), avg_stars, stars.size());

        bool quit = false;
        if (pid == 0) {
            quit = check_quit();
        }
        
        if (use_mmpi) {
            mmpi_bcast(0, &quit, sizeof(bool)); // TODO: this is probably also adding some latency
        } else {
            MPI_Bcast(&quit, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
        }
        if (quit) break;
    }
    auto run_end = chrono::now();
    millis run_time = run_end - run_start;
    if (pid == 0) display_cleanup();

    if (use_mmpi) {
        mmpi_finalize();
    } else {
        MPI_Finalize();
    }

    if (pid == 0) {
        fprintf(stdout, "%d iterations took %.01fms for %.01fms each\n", 
            NUM_ITERS, run_time.count(), run_time.count()/NUM_ITERS);
    }
}

/**
 * @brief Entry point for distributed mpi implementation
 */
void mpi_distributed_main(int argc, char* argv[]) {

    use_mmpi = true;

    printf("Hello MPI distributed\n");

    if (argc != 4) {
        fprintf(stderr, "Need 3 arguments for distributed\n");
        exit(-1);
    }

    int init_pid = atoi(argv[2]);
    int init_node_ct = atoi(argv[3]);

    mmpi_init(init_pid, init_node_ct);
    pid = init_pid;
    nprocs = init_node_ct;

    mpi_run_simulation();
}

/**
 * @brief Entry point for single-node MPI implementation
 */
void mpi_single_main(int argc, char* argv[]) {

    use_mmpi = false;

    fprintf(stderr, "Hello MPI single\n");

#ifndef USE_MPI

    fprintf(stderr, "MPI not enabled!\n");
    exit(-1);

#endif /* USE_MPI */

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    mpi_run_simulation();
}