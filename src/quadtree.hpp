/**
 * @file quadtree.hpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Header file for quadtree implementation
 */

#pragma once

#include "compact_defines.h"
#include "stars.hpp"
#include <vector>

struct QNode {
    float cx, cy;     // center of region
    float side_len;

    QNode* nw = nullptr;
    QNode* ne = nullptr;
    QNode* sw = nullptr;
    QNode* se = nullptr;

    Star* s = nullptr;
    float mass = 0.0f;
    float com_x = 0.0f; // center of mass
    float com_y = 0.0f;
};

QNode* build_qtree(std::vector<Star> &stars);
void destroy_tree(QNode* node);