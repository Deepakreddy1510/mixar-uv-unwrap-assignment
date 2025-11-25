/**
 * @file topology.h
 * @brief Mesh topology (edges and adjacency)
 *
 * PROVIDED API - YOU IMPLEMENT IN topology.cpp
 */

#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Topology information for a mesh
 *
 * Stores:
 * - All unique edges in the mesh
 * - For each edge, the 1 or 2 adjacent faces
 */
typedef struct {
    int* edges;            /**< Edge vertex pairs [v0,v1, v0,v1, ...] (2 * num_edges) */
    int num_edges;         /**< Number of edges */

    int* edge_faces;       /**< Adjacent faces [f0,f1, f0,f1, ...] (2 * num_edges)
                                 f1 = -1 for boundary edges */
} TopologyInfo;

/**
 * @brief Build topology from mesh
 *
 * Algorithm:
 * 1. Iterate through all triangles
 * 2. Extract all edges, ensuring uniqueness (v0 < v1)
 * 3. For each edge, find adjacent faces (1 or 2)
 * 4. Validate using Euler characteristic
 *
 * @param mesh Input mesh
 * @return Newly allocated topology info, or NULL on error
 * @note Caller must free with free_topology()
 *
 * IMPLEMENTATION REQUIRED
 */
TopologyInfo* build_topology(const Mesh* mesh);

/**
 * @brief Free topology memory
 * @param topo Topology to free
 */
void free_topology(TopologyInfo* topo);

/**
 * @brief Validate topology using Euler characteristic
 * @param mesh Original mesh
 * @param topo Computed topology
 * @return 1 if valid, 0 otherwise
 *
 * For closed mesh: V - E + F = 2
 */
int validate_topology(const Mesh* mesh, const TopologyInfo* topo);

#ifdef __cplusplus
}
#endif

#endif /* TOPOLOGY_H */
