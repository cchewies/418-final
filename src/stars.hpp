#pragma once
#include "compact_defines.h"
#include <vector>

struct Star {
    float x, y;     // position
    float vx, vy;   // velocity
    float mass;     // mass
};

std::vector<Star> generate_galaxy(void);
