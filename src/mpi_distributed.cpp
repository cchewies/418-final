#include "mpi_distributed.hpp"
#include "compact_defines.h"
#include "minimpi.hpp"
#include "simulation_config.hpp"
#include "display.hpp"
#include "quadtree.hpp"
#include <cmath>
#include <algorithm>

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

    int my_start, my_end;

    // // compute spatial partition (pure naive)
    // my_start = stars.size() * mmpi_getpid() / mmpi_getnodes();
    // my_end   = stars.size() * (mmpi_getpid() + 1) / mmpi_getnodes();

    if (stars[0].cost == 0) { // first iter
        // naive allocation
        my_start = stars.size() * (mmpi_getpid()) / mmpi_getnodes();
        my_end = stars.size() * (mmpi_getpid()+1) / mmpi_getnodes();
    } else {
        // prefix-sum load balancing

        std::vector<int> prefix(stars.size());
        prefix[0] = stars[0].cost;
        for (size_t i = 1; i < stars.size(); i++) {
            prefix[i] = prefix[i-1] + stars[i].cost;
        }

        int total = prefix.back();

        // compute boundaries up front for consistency 
        std::vector<int> boundaries(mmpi_getnodes() + 1);
        boundaries[0] = 0;
        boundaries[mmpi_getnodes()] = stars.size();

        for (int r = 1; r < mmpi_getnodes(); r++) {
            int target = (int)((float)r / mmpi_getnodes() * total);
            // Find first index where prefix >= target
            boundaries[r] = (int)(std::lower_bound(prefix.begin(), prefix.end(), target) 
                                - prefix.begin());
        }

        my_start = boundaries[mmpi_getpid()];
        my_end   = boundaries[mmpi_getpid() + 1];
    }

    // -- end of stars assignment --

    auto qtree_start = chrono::now();
    QNode* root = build_qtree(stars);
    auto qtree_end = chrono::now();
    millis qtree_time = qtree_end - qtree_start;

    int total_ct = 0;

    // // get node's allocation
    // int my_start = NUM_STARS * (mmpi_getpid()) / mmpi_getnodes();
    // int my_end = NUM_STARS * (mmpi_getpid()+1) / mmpi_getnodes();

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
    std::vector<int> counts(mmpi_getnodes()), displs(mmpi_getnodes());

    int my_count = (my_end - my_start) * sizeof(Star);
    // MPI_Allgather(&my_count, 1, MPI_INT, counts.data(), 1, MPI_INT, MPI_COMM_WORLD); // TODO: convert
    counts[mmpi_getpid()] = my_count;
    for (int r = 0; r < mmpi_getnodes(); r++) {
        mmpi_bcast(r, &counts[r], sizeof(int));
    }

    // for (int vpid = 0; vpid < mmpi_getnodes(); vpid++) {
    //     int node_start = NUM_STARS * vpid / mmpi_getnodes();
    //     int node_end = NUM_STARS * (vpid + 1) / mmpi_getnodes();
    //     counts[vpid] = (node_end - node_start) * sizeof(Star);
    // }

    displs[0] = 0;
    for (int vpid = 1; vpid < mmpi_getnodes(); vpid++) {
        displs[vpid] = displs[vpid - 1] + counts[vpid - 1];
    }

    // gather updated stars to all ranks
    // std::vector<Star> recv_buf(stars.size());
    // mmpi_syncv(recv_buf.data(), stars.size() * sizeof(Star), counts.data(), displs.data()); // TODO: check
    // stars = std::move(recv_buf);
    // test -start-
    std::vector<Star> recv_buf(stars.size());
    std::memcpy(
        recv_buf.data() + my_start,
        stars.data() + my_start,
        (my_end - my_start) * sizeof(Star)
    );
    mmpi_syncv(
        recv_buf.data(),
        stars.size() * sizeof(Star),
        counts.data(),
        displs.data()
    );
    stars = std::move(recv_buf);
    // test -end-

    auto comm_end = chrono::now();
    millis comm_time = comm_end - comm_start;

    fprintf(stdout, "Qtree took %.01fms, force took %.01fms, comm took %.01fms\n", 
            qtree_time.count(), force_time.count(), comm_time.count());

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
