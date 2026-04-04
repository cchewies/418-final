/**
 * @file display.hpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Display config
 */

#pragma once

#include "compact_defines.h"
#include "stars.hpp"

#define RENDER_ENABLED

// Window size
constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;

// rendering period (because X11 forwarding is slow)
constexpr int RENDER_PERIOD = 1000;

void display_init(void);
void display_check_key(void);
void display_render(std::vector<Star> &stars);