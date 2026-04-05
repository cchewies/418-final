/**
 * @file stars.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Various initial starting configurations
 */

#include "stars.hpp"
#include "display.hpp"
#include "simulation_config.hpp"
#include <vector>
#include <ctime>
#include <math.h>

/**
 * @brief Random blob galaxy
 */
static std::vector<Star> random_blob(void) {
    std::vector<Star> stars;
    std::srand(static_cast<unsigned>(std::time(nullptr)));

    // gravitational constant
    const float center_x = WINDOW_WIDTH / 2.0f;
    const float center_y = WINDOW_HEIGHT / 2.0f;
    const float galaxy_radius = std::min(WINDOW_WIDTH, WINDOW_HEIGHT) / 2.5f;

    // center black hole
    Star black_hole;
    black_hole.x = center_x;
    black_hole.y = center_y;
    black_hole.vx = 0;
    black_hole.vy = 0;
    black_hole.mass = 1000.0f;
    stars.push_back(black_hole);

    for (int i = 0; i < NUM_STARS - 1; i++) {
        Star s;

        // Random radius with denser core
        float r = galaxy_radius * (float)std::rand() / RAND_MAX;
        float angle = 2.0f * M_PI * ((float)std::rand() / RAND_MAX);

        // Polar to Cartesian
        s.x = center_x + r * cos(angle);
        s.y = center_y + r * sin(angle);

        // random mass between 1 and 2
        s.mass = 1.0f + ((float)(std::rand() % 10) / 10.0f);

        // tangent velocity + random bits
        float dx = s.x - center_x;
        float dy = s.y - center_y;
        float dist = sqrt(dx*dx + dy*dy);
        float v_circ = (1+3*(r*r/galaxy_radius/galaxy_radius)) * sqrt(G * black_hole.mass / (dist + 1.0f));
        s.vx = -dy / dist * v_circ;
        s.vy = dx / dist * v_circ;
        s.vx += ((std::rand() % 20) - 10) / 200.0f;
        s.vy += ((std::rand() % 20) - 10) / 200.0f;

        stars.push_back(s);
    }
    return stars;
}

/**
 * @brief Generate a random galaxy
 * 
 * @return Vector of stars
 */
std::vector<Star> generate_galaxy(void) {
    return random_blob();
}
