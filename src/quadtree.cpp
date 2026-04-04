#include "quadtree.hpp"

#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>

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
    if (s->x < node->cx) {
        if (s->y < node->cy) return node->sw;
        else return node->nw;
    } else {
        if (s->y < node->cy) return node->se;
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
    node->com_x = (node->com_x * node->mass + s->x * s->mass) / new_mass;
    node->com_y = (node->com_y * node->mass + s->y * s->mass) / new_mass;
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

    // compute bounding box
    float min_x = stars[0].x, max_x = stars[0].x;
    float min_y = stars[0].y, max_y = stars[0].y;

    for (auto& s : stars) {
        min_x = std::min(min_x, s.x);
        max_x = std::max(max_x, s.x);
        min_y = std::min(min_y, s.y);
        max_y = std::max(max_y, s.y);
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
}