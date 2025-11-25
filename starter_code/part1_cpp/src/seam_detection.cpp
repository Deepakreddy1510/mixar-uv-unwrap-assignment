/**
 * @file seam_detection.cpp
 * @brief Seam detection using spanning tree + angular defect
 *
 * SKELETON - YOU IMPLEMENT THIS
 *
 * Algorithm:
 * 1. Build dual graph (faces as nodes, shared edges as edges)
 * 2. Compute spanning tree via BFS
 * 3. Mark non-tree edges as seam candidates
 * 4. Refine using angular defect
 *
 * See reference/algorithms.md for detailed description
 */

#include "unwrap.h"
#include "math_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <set>
#include <queue>

/**
 * @brief Compute angular defect at a vertex
 *
 * Angular defect = 2π - sum of angles at vertex
 *
 * - Flat surface: defect ≈ 0
 * - Corner (like cube): defect > 0
 * - Saddle: defect < 0
 *
 * @param mesh Input mesh
 * @param vertex_idx Vertex index
 * @return Angular defect in radians
 */
static float compute_angular_defect(const Mesh* mesh, int vertex_idx) {
    // TODO: Implement
    //
    // Steps:
    // 1. Find all triangles containing this vertex
    // 2. For each triangle, compute angle at this vertex
    //    (use compute_vertex_angle_in_triangle from math_utils.h)
    // 3. Sum all angles
    // 4. Return 2*PI - sum
    //
    // Hint: Angular defect indicates curvature
    //       High defect → sharp feature → good seam location

    float angle_sum = 0.0f;

    // YOUR CODE HERE

    return 2.0f * M_PI - angle_sum;
}

/**
 * @brief Get all edges incident to a vertex
 */
static std::vector<int> get_vertex_edges(const TopologyInfo* topo, int vertex_idx) {
    std::vector<int> edges;

    // TODO: Iterate through all edges, add those touching vertex_idx

    // YOUR CODE HERE

    return edges;
}

int* detect_seams(const Mesh* mesh,
                  const TopologyInfo* topo,
                  float angle_threshold,
                  int* num_seams_out) {
    if (!mesh || !topo || !num_seams_out) return NULL;

    // TODO: Implement seam detection
    //
    // ALGORITHM:
    //
    // STEP 1: Build dual graph
    //   - Nodes = faces
    //   - Edges = shared edges between faces
    //   - Use std::vector<std::vector<int>> for adjacency list
    //
    // STEP 2: Spanning tree via BFS
    //   - Start from face 0
    //   - Mark edges in spanning tree
    //   - Use std::set<int> to track tree edges
    //
    // STEP 3: Initial seam candidates = non-tree edges
    //   - All edges NOT in spanning tree
    //
    // STEP 4: Angular defect refinement
    //   - For each vertex with high angular defect (> 0.5 radians)
    //   - Add incident edges to seam candidates
    //
    // STEP 5: Convert seam candidates to array
    //
    // Expected seam counts:
    //   Cube: 7-9 seams
    //   Sphere: 1-3 seams
    //   Cylinder: 1-2 seams

    std::set<int> seam_candidates;

    // YOUR CODE HERE

    // Convert to array
    *num_seams_out = seam_candidates.size();
    int* seams = (int*)malloc(*num_seams_out * sizeof(int));

    int idx = 0;
    for (int edge_idx : seam_candidates) {
        seams[idx++] = edge_idx;
    }

    printf("Detected %d seams\n", *num_seams_out);

    return seams;
}
