#include "compact_defines.h"
#include "display.hpp"
#include "serial.hpp"
#include "simulation_config.hpp"
#include <math.h>
#include <chrono>

/**
 * @brief Naive N^2 serial implementation
 * 
 * @param stars 
 */
static void iterate_simulation(std::vector<Star> &stars) {

    // calculate pairwise forces
    for (size_t i = 0; i < stars.size(); i++) {
        float fx = 0, fy = 0;
        for (size_t j = 0; j < stars.size(); j++) {
            if (i == j) continue;
            float dx = stars[j].x - stars[i].x;
            float dy = stars[j].y - stars[i].y;
            float dist2 = dx*dx + dy*dy + 1;
            float dist = sqrt(dist2);
            float force = G * stars[i].mass * stars[j].mass / dist2;
            fx += force * dx / dist;
            fy += force * dy / dist;
        }
        stars[i].vx += fx / stars[i].mass * dt;
        stars[i].vy += fy / stars[i].mass * dt;
    }

    // update positions
    for (auto& s : stars) {
        s.x += s.vx * dt;
        s.y += s.vy * dt;
    }
}

/**
 * @brief Entry point for serial implementation
 */
void serial_main(void) {

    std::vector<Star> stars = generate_galaxy();

    float slowest_frame = 0;
    while (1) {
        for (int i = 0; i < 100; i++) {
            auto start = std::chrono::high_resolution_clock::now();
    
            iterate_simulation(stars);
            
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::milli> frame_time = end - start;

            if (frame_time.count() > slowest_frame) {
                slowest_frame = frame_time.count();
            }
        }
        display_render(stars);
        fprintf(stdout, "Slowest iteration: %f\n", slowest_frame);
        slowest_frame = 0;
    }
}
