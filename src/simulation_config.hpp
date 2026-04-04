#pragma once

// star count (+1 for black hole)
constexpr int NUM_STARS = 100000;

// Simulation time interval
constexpr float DT = 0.001;

// Gravitational constant
constexpr float G = 10;

// Barnes–Hut threshold
constexpr float THETA = 0.5f;

// softening
constexpr float EPS2 = 1;

// rendering period (because X11 forwarding is slow)
constexpr int RENDER_PERIOD = 1000;