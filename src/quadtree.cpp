/**
 * @file quadtree.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Quadtree implementation for stars
 */

#include "quadtree.hpp"
#include "simulation_config.hpp"

#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>

static int num_alloc = 0;

/**
 * @brief Create a quadtree node 
 * 
 * @param cx Center x coordinate
 * @param cy Center y coordinate
 * @param side_len Side length of cell
 * @return Created node
 */
static QNode* create_node(float cx, float cy, float side_len) {
    QNode* node = new QNode{};
    num_alloc++;
    node->cx = cx;
    node->cy = cy;
    node->side_len = side_len;
    return node;
}

/**
 * @brief Subdivide a node into 4 quadrants
 * 
 * @param node Node to subdivide
 */
static void subdivide(QNode* node) {
    float half = node->side_len / 2.0f;

    node->ne = create_node(node->cx + half/2, node->cy + half/2, half);
    node->sw = create_node(node->cx - half/2, node->cy - half/2, half);
    node->nw = create_node(node->cx - half/2, node->cy + half/2, half);
    node->se = create_node(node->cx + half/2, node->cy - half/2, half);
}

/**
 * @brief Get a reference to the child a star would belong to
 * 
 * @param node Node to get child of
 * @param s Star to get correct child for
 * @return Reference to the child cell
 */
static QNode*& get_child(QNode* node, Star* s) {
    if (s->pos.x < node->cx) {
        if (s->pos.y < node->cy) return node->sw;
        else return node->nw;
    } else {
        if (s->pos.y < node->cy) return node->se;
        else return node->ne;
    }
}

/**
 * @brief Insert a star at this node (may subdivide)
 * 
 * @param node Node to insert into
 * @param s Star to insert
 */
static void insert_star(QNode* node, Star* s) {
    // update mass and center of mass
    float new_mass = node->mass + s->mass;
    assert(new_mass != 0);  // we dont have 0 mass stars
    node->com_x = (node->com_x * node->mass + s->pos.x * s->mass) / new_mass;
    node->com_y = (node->com_y * node->mass + s->pos.y * s->mass) / new_mass;
    node->mass = new_mass;

    // empty leaf, just insert star
    if (node->s == nullptr && node->nw == nullptr) {
        node->s = s;
        return;
    }

    // leaf but already has a star -> subdivide
    if (node->s != nullptr && node->nw == nullptr) {
        subdivide(node);

        // reinsert existing star
        Star* old = node->s;
        node->s = nullptr;
        insert_star(get_child(node, old), old);
    }

    // insert new star
    insert_star(get_child(node, s), s);
}

/**
 * @brief Build a quadtree from a list of stars
 * 
 * @param stars Stars list (should not be empty)
 * @return Constructed quadtree
 */
QNode* build_qtree(std::vector<Star>& stars) {
    assert(!stars.empty());
    assert(num_alloc == 0);

    // compute bounding box
    float min_x = stars[0].pos.x, max_x = stars[0].pos.x;
    float min_y = stars[0].pos.y, max_y = stars[0].pos.y;

    for (auto& s : stars) {
        min_x = std::min(min_x, s.pos.x);
        max_x = std::max(max_x, s.pos.x);
        min_y = std::min(min_y, s.pos.y);
        max_y = std::max(max_y, s.pos.y);
    }

    float cx = (min_x + max_x) / 2.0f;
    float cy = (min_y + max_y) / 2.0f;
    float side_len = std::max(max_x - min_x, max_y - min_y);

    // padding to avoid boundary issues
    side_len += 1e-5f;

    QNode* root = create_node(cx, cy, side_len);

    // add all stars
    for (auto& s : stars) {
        insert_star(root, &s);
    }

    return root;
}

/**
 * @brief Free resources used by a quadtree node
 * 
 * @param node Node to free
 */
void destroy_tree(QNode* node) {
    if (!node) return;

    destroy_tree(node->nw);
    destroy_tree(node->ne);
    destroy_tree(node->sw);
    destroy_tree(node->se);

    delete node;
    num_alloc--;
}

// expand 32-bit int into 64-bit w interleaved zeros
u64 expand_bits(u32 x) {
    u64 v = x;
    v = (v | (v << 16)) & 0x0000FFFF0000FFFFULL;
    v = (v | (v << 8))  & 0x00FF00FF00FF00FFULL;
    v = (v | (v << 4))  & 0x0F0F0F0F0F0F0F0FULL;
    v = (v | (v << 2))  & 0x3333333333333333ULL;
    v = (v | (v << 1))  & 0x5555555555555555ULL;
    return v;
}

// 2D Morton encoding
u64 morton2D(float x, float y, float min_x, float min_y, float max_x, float max_y) {
    float nx = (x - min_x) / (max_x - min_x); // normalize to [0, 1]
    float ny = (y - min_y) / (max_y - min_y);
    nx = std::min(1.0f, std::max(0.0f, nx));
    ny = std::min(1.0f, std::max(0.0f, ny));
    u32 ix = (u32)(nx * ((1 << 21) - 1));
    u32 iy = (u32)(ny * ((1 << 21) - 1));
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
int compute_force(Star& s, QNode* node, float& fx, float& fy) {

    int star_count = 0;

    // Node doesnt exist, no stars, no children, or star is current one
    if (!node || node->mass == 0 || (node->s == &s && !node->nw)) {
        return star_count;
    }

    float dx = node->com_x - s.pos.x;
    float dy = node->com_y - s.pos.y;
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