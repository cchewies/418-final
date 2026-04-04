/**
 * @file serial.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Single-threaded Barnes-Hut implementations
 * 
 * naive: O(N^2), single-thread, single-node
 * serial: O(NlogN), single-thread, single-node
 */

#include "compact_defines.h"
#include "display.hpp"
#include "serial.hpp"
#include "simulation_config.hpp"
#include "quadtree.hpp"
#include <math.h>
#include <chrono>

/**
 * @brief Naive N^2 serial implementation
 * 
 * @param stars 
 */
static void naive_iterate_simulation(std::vector<Star> &stars) {

    // calculate pairwise forces
    for (size_t i = 0; i < stars.size(); i++) {
        float fx = 0, fy = 0;
        for (size_t j = 0; j < stars.size(); j++) {
            if (i == j) continue;
            float dx = stars[j].x - stars[i].x;
            float dy = stars[j].y - stars[i].y;
            float dist2 = dx*dx + dy*dy + EPS2;
            float dist = sqrt(dist2);
            float force = G * stars[i].mass * stars[j].mass / dist2;
            fx += force * dx / dist;
            fy += force * dy / dist;
        }
        stars[i].vx += fx / stars[i].mass * DT;
        stars[i].vy += fy / stars[i].mass * DT;
    }

    // update positions
    for (auto& s : stars) {
        s.x += s.vx * DT;
        s.y += s.vy * DT;
    }
}

/**
 * @brief Entry point for naive implementation
 */
void naive_main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;

    fprintf(stderr, "Hello naive\n");

    std::vector<Star> stars = generate_galaxy();

    while (1) {
        display_render(stars);

        auto start = std::chrono::high_resolution_clock::now();

        naive_iterate_simulation(stars);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> frame_time = end - start;

        fprintf(stdout, "Iteration took %.00fms, stars: %ld\n", 
            frame_time.count(), stars.size());
        display_check_key();
    }
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
int compute_force(Star& s, QNode* node, float& fx, float& fy) {

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
static int serial_iterate_simulation(std::vector<Star> &stars) {

    QNode* root = build_qtree(stars);

    int total_ct = 0;

    // update velocities
    for (auto& s : stars) {
        float fx = 0.0f;
        float fy = 0.0f;

        int star_ct = compute_force(s, root, fx, fy);
        total_ct += star_ct;

        s.vx += fx / s.mass * DT;
        s.vy += fy / s.mass * DT;
    }

    // update positions
    for (auto& s : stars) {
        s.x += s.vx * DT;
        s.y += s.vy * DT;
    }

    destroy_tree(root);
    return (int)(total_ct / stars.size());
}

/**
 * @brief Entry point for serial quadtree implementation
 */
void serial_main(int argc, char* argv[]) {
    (void)argc;
    (void)argv;

    fprintf(stderr, "Hello serial\n");

    std::vector<Star> stars = generate_galaxy();

    while (1) {
        display_render(stars);
        auto start = std::chrono::high_resolution_clock::now();

        int avg_stars = serial_iterate_simulation(stars);
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> frame_time = end - start;

        fprintf(stdout, "Iteration took %.00fms, avg stars: %d/%ld\n", 
            frame_time.count(), avg_stars, stars.size());
        display_check_key();
    }
}
