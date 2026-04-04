#pragma once

#include "compact_defines.h"
#include "stars.hpp"

// Window size
constexpr int WINDOW_WIDTH = 1280;
constexpr int WINDOW_HEIGHT = 720;

void display_init(void);
void display_check_key(void);
void display_render(std::vector<Star> &stars);