/**
 * @file simulation_config.hpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Config for simulation
 */

#pragma once

// star count
constexpr int NUM_STARS = 5000;

// Simulation time interval
constexpr float DT = 0.001;

// Gravitational constant
constexpr float G = 10;

// Barnes-Hut threshold
constexpr float THETA = 0.5f;

// softening
constexpr float EPS2 = 1;
