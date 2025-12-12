/**
 * @file topology.cpp
 * @brief Topology builder implementation
 *
 * SKELETON - YOU IMPLEMENT THIS
 *
 * Algorithm:
 * 1. Extract all edges from triangles
 * 2. Ensure uniqueness (always store as v0 < v1)
 * 3. For each edge, find adjacent faces
 * 4. Validate using Euler characteristic
 */

#include "topology.h"
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <vector>

/**
 * @brief Edge structure for uniqueness
 */
struct Edge {
    int v0, v1;

    Edge(int a, int b) {
        // Always store smaller vertex first
        if (a < b) {
            v0 = a;
            v1 = b;
        } else {
            v0 = b;
            v1 = a;
        }
    }

    bool operator<(const Edge& other) const {
        if (v0 != other.v0) return v0 < other.v0;
        return v1 < other.v1;
    }
};

/**
 * @brief Edge information
 */
struct EdgeInfo {
    int face0;
    int face1;

    EdgeInfo() : face0(-1), face1(-1) {}
};

TopologyInfo* build_topology(const Mesh* mesh) {
    if (!mesh) return NULL;

    // TODO: Implement topology building
    //
    // Steps:
    // 1. Create std::map<Edge, EdgeInfo> to collect edges
    // 2. Iterate through all triangles
    //    For each triangle, extract 3 edges
    //    Add to map, tracking which faces use each edge
    // 3. Convert map to arrays (edges, edge_faces)
    // 4. Allocate TopologyInfo and fill arrays
    //
    // Hints:
    // - Use Edge struct for automatic ordering
    // - Each edge should have 1 or 2 adjacent faces
    // - Boundary edges have only 1 face (set face1 = -1)
    //
    // See reference/topology_example.cpp for complete example

    TopologyInfo* topo = (TopologyInfo*)malloc(sizeof(TopologyInfo));

    // Initialize to safe defaults (prevents crashes before implementation)
    topo->edges = NULL;
    topo->num_edges = 0;
    topo->edge_faces = NULL;

    // ===============================
    // BUILD TOPOLOGY IMPLEMENTATION
    // ===============================
    std::map<Edge, EdgeInfo> edge_map;

    // 1. Collect edges from all triangles (mesh->triangles is flat: 3 * num_triangles)
    for (int f = 0; f < mesh->num_triangles; ++f) {
        int base = 3 * f;
        int v0 = mesh->triangles[base + 0];
        int v1 = mesh->triangles[base + 1];
        int v2 = mesh->triangles[base + 2];

        Edge e0(v0, v1);
        Edge e1(v1, v2);
        Edge e2(v2, v0);

        Edge edgesArr[3] = { e0, e1, e2 };
        for (int i = 0; i < 3; ++i) {
            Edge &Ekey = edgesArr[i];
            EdgeInfo &info = edge_map[Ekey]; // default-constructed if missing
            if (info.face0 == -1) info.face0 = f;
            else if (info.face1 == -1) info.face1 = f;
            // If there are more than 2 adjacent faces something is wrong; ignore extra.
        }
    }

    // 2. Convert map to plain arrays (TopologyInfo expects int* arrays)
    int numE = (int)edge_map.size();
    topo->num_edges = numE;

    if (numE > 0) {
        topo->edges = (int*)malloc(sizeof(int) * 2 * numE);
        topo->edge_faces = (int*)malloc(sizeof(int) * 2 * numE);
        if (!topo->edges || !topo->edge_faces) {
            // allocation failed - cleanup and return NULL
            if (topo->edges) free(topo->edges);
            if (topo->edge_faces) free(topo->edge_faces);
            free(topo);
            return NULL;
        }
    } else {
        topo->edges = NULL;
        topo->edge_faces = NULL;
    }

    int idx = 0;
    for (const auto &kv : edge_map) {
        const Edge &e = kv.first;
        const EdgeInfo &info = kv.second;

        topo->edges[2*idx + 0] = e.v0;
        topo->edges[2*idx + 1] = e.v1;

        topo->edge_faces[2*idx + 0] = info.face0;
        topo->edge_faces[2*idx + 1] = info.face1;

        ++idx;
    }



    return topo;
}

void free_topology(TopologyInfo* topo) {
    if (!topo) return;

    if (topo->edges) free(topo->edges);
    if (topo->edge_faces) free(topo->edge_faces);
    free(topo);
}

int validate_topology(const Mesh* mesh, const TopologyInfo* topo) {
    if (!mesh || !topo) return 0;

    int V = mesh->num_vertices;
    int E = topo->num_edges;
    int F = mesh->num_triangles;

    int euler = V - E + F;

    printf("Topology validation:\n");
    printf("  V=%d, E=%d, F=%d\n", V, E, F);
    printf("  Euler characteristic: %d (expected 2 for closed mesh)\n", euler);

    // Closed meshes should have Euler = 2
    // Open meshes or meshes with holes may differ
    if (euler != 2) {
        printf("  Warning: Non-standard Euler characteristic\n");
        printf("  (This may be OK for open meshes or meshes with boundaries)\n");
    }

    return 1;
}
