/**
 * @file unwrap.cpp
 * @brief Main UV unwrapping orchestrator
 *
 * SKELETON - YOU IMPLEMENT THIS
 *
 * This file ties together all the components:
 * - Topology building
 * - Seam detection
 * - Island extraction
 * - LSCM parameterization
 * - Island packing
 */

#include "unwrap.h"
#include "lscm.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <map>

/**
 * @brief Extract UV islands after seam cuts
 *
 * Uses connected components algorithm on face graph after removing seam edges.
 *
 * @param mesh Input mesh
 * @param topo Topology info
 * @param seam_edges Array of seam edge indices
 * @param num_seams Number of seams
 * @param num_islands_out Output: number of islands
 * @return Array of island IDs per face
 */
static int* extract_islands(const Mesh* mesh,
                           const TopologyInfo* topo,
                           const int* seam_edges,
                           int num_seams,
                           int* num_islands_out) {
    // TODO: Implement island extraction
    //
    // Algorithm:
    // 1. Build face adjacency graph (faces connected by non-seam edges)
    // 2. Run connected components (BFS or DFS)
    // 3. Assign island ID to each face
    //
    // Steps:
    // 1. Create std::set<int> of seam edge indices for fast lookup
    // 2. Build adjacency list for faces (only through non-seam edges)
    // 3. Run BFS/DFS to find connected components
    // 4. Return array of island IDs (one per face)

    int* face_island_ids = (int*)malloc(mesh->num_triangles * sizeof(int));

    // Initialize all to -1 (unvisited)
    for (int i = 0; i < mesh->num_triangles; i++) {
        face_island_ids[i] = -1;
    }

    // YOUR CODE HERE

    *num_islands_out = 0;  // Update this with actual count

    printf("Extracted %d UV islands\n", *num_islands_out);

    return face_island_ids;
}

/**
 * @brief Copy UVs from island parameterization to result mesh
 */
static void copy_island_uvs(Mesh* result,
                           const float* island_uvs,
                           const int* face_indices,
                           int num_faces,
                           const std::map<int, int>& global_to_local) {
    // TODO: Implement UV copying
    //
    // For each face in the island:
    //   For each vertex in the face:
    //     global_idx = vertex index in mesh
    //     local_idx = global_to_local[global_idx]
    //     result->uvs[global_idx * 2] = island_uvs[local_idx * 2]
    //     result->uvs[global_idx * 2 + 1] = island_uvs[local_idx * 2 + 1]

    // YOUR CODE HERE
}

Mesh* unwrap_mesh(const Mesh* mesh,
                  const UnwrapParams* params,
                  UnwrapResult** result_out) {
    if (!mesh || !params || !result_out) {
        fprintf(stderr, "unwrap_mesh: Invalid arguments\n");
        return NULL;
    }

    printf("\n=== UV Unwrapping ===\n");
    printf("Input: %d vertices, %d triangles\n",
           mesh->num_vertices, mesh->num_triangles);
    printf("Parameters:\n");
    printf("  Angle threshold: %.1f°\n", params->angle_threshold);
    printf("  Min island faces: %d\n", params->min_island_faces);
    printf("  Pack islands: %s\n", params->pack_islands ? "yes" : "no");
    printf("  Island margin: %.3f\n", params->island_margin);
    printf("\n");

    // TODO: Implement main unwrapping pipeline
    //
    // STEP 1: Build topology
    TopologyInfo* topo = build_topology(mesh);
    if (!topo) {
        fprintf(stderr, "Failed to build topology\n");
        return NULL;
    }
    validate_topology(mesh, topo);

    // STEP 2: Detect seams
    int num_seams;
    int* seam_edges = detect_seams(mesh, topo, params->angle_threshold, &num_seams);

    // STEP 3: Extract islands
    int num_islands;
    int* face_island_ids = extract_islands(mesh, topo, seam_edges, num_seams, &num_islands);

    // STEP 4: Parameterize each island using LSCM
    Mesh* result = allocate_mesh_copy(mesh);
    result->uvs = (float*)calloc(mesh->num_vertices * 2, sizeof(float));

    for (int island_id = 0; island_id < num_islands; island_id++) {
        printf("\nProcessing island %d/%d...\n", island_id + 1, num_islands);

        // Get faces in this island
        std::vector<int> island_faces;
        for (int f = 0; f < mesh->num_triangles; f++) {
            if (face_island_ids[f] == island_id) {
                island_faces.push_back(f);
            }
        }

        printf("  %d faces in island\n", (int)island_faces.size());

        if (island_faces.size() < params->min_island_faces) {
            printf("  Skipping (too small)\n");
            continue;
        }

        // YOUR CODE HERE:
        // - Call lscm_parameterize
        // - Build global_to_local mapping
        // - Copy UVs to result mesh
    }

    // STEP 5: Pack islands if requested
    if (params->pack_islands) {
        UnwrapResult temp_result;
        temp_result.num_islands = num_islands;
        temp_result.face_island_ids = face_island_ids;

        pack_uv_islands(result, &temp_result, params->island_margin);
    }

    // STEP 6: Compute quality metrics
    UnwrapResult* result_data = (UnwrapResult*)malloc(sizeof(UnwrapResult));
    result_data->num_islands = num_islands;
    result_data->face_island_ids = face_island_ids;
    compute_quality_metrics(result, result_data);

    *result_out = result_data;

    // Cleanup
    free_topology(topo);
    free(seam_edges);

    printf("\n=== Unwrapping Complete ===\n");

    return result;
}

void free_unwrap_result(UnwrapResult* result) {
    if (!result) return;

    if (result->face_island_ids) {
        free(result->face_island_ids);
    }
    free(result);
}
