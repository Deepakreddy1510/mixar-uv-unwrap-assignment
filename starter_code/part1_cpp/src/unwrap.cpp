/**
 * @file unwrap.cpp
 * @brief Main UV unwrapping orchestrator
 *
 * IMPLEMENTATION PROVIDED BELOW
 */

#include "unwrap.h"
#include "lscm.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <string.h> // for memcpy

/**
 * @brief Extract UV islands after seam cuts
 *
 * Uses connected components algorithm on face graph after removing seam edges.
 */
static int* extract_islands(const Mesh* mesh,
                           const TopologyInfo* topo,
                           const int* seam_edges,
                           int num_seams,
                           int* num_islands_out) {
    if (!mesh || !topo || !num_islands_out) return NULL;

    int F = mesh->num_triangles;
    int E = topo->num_edges;

    // Build set of seam edges for quick lookup
    std::set<int> seam_set;
    for (int i = 0; i < num_seams; ++i) {
        seam_set.insert(seam_edges[i]);
    }

    // Build face adjacency list (without crossing seam edges)
    std::vector<std::vector<int>> adj(F);
    for (int ei = 0; ei < E; ++ei) {
        int f0 = topo->edge_faces[2*ei + 0];
        int f1 = topo->edge_faces[2*ei + 1];

        // Only consider internal edges with two adjacent faces
        if (f0 >= 0 && f1 >= 0) {
            // If this edge is a seam, skip connecting the two faces
            if (seam_set.find(ei) != seam_set.end()) continue;

            adj[f0].push_back(f1);
            adj[f1].push_back(f0);
        }
    }

    // Prepare output array (initialized to -1)
    int* face_island_ids = (int*)malloc(sizeof(int) * F);
    if (!face_island_ids) return NULL;
    for (int i = 0; i < F; ++i) face_island_ids[i] = -1;

    int island_id = 0;
    std::vector<int> queue_buf;
    for (int f = 0; f < F; ++f) {
        if (face_island_ids[f] != -1) continue; // already assigned

        // BFS/DFS from face f
        queue_buf.clear();
        queue_buf.push_back(f);
        face_island_ids[f] = island_id;

        size_t qhead = 0;
        while (qhead < queue_buf.size()) {
            int cur = queue_buf[qhead++];
            for (int nb : adj[cur]) {
                if (face_island_ids[nb] == -1) {
                    face_island_ids[nb] = island_id;
                    queue_buf.push_back(nb);
                }
            }
        }

        island_id++;
    }

    *num_islands_out = island_id;
    printf("Extracted %d UV islands\n", *num_islands_out);
    return face_island_ids;
}

/**
 * @brief Copy UVs from island parameterization to result mesh
 *
 * - island_uvs: pointer to 2 * num_local_vertices floats [u0,v0,u1,v1,...]
 * - face_indices: array of global face indices that belong to the island
 * - global_to_local: mapping from global vertex index -> local index in island_uvs
 */
static void copy_island_uvs(Mesh* result,
                           const float* island_uvs,
                           const int* face_indices,
                           int num_faces,
                           const std::map<int, int>& global_to_local) {
    if (!result || !island_uvs || !face_indices) return;

    // For each face, write UVs for each of its 3 vertices
    for (int fi = 0; fi < num_faces; ++fi) {
        int face_idx = face_indices[fi];
        int base = 3 * face_idx;
        int g0 = result->triangles[base + 0];
        int g1 = result->triangles[base + 1];
        int g2 = result->triangles[base + 2];

        // For each global vertex, if we have a mapping, copy UV
        auto it0 = global_to_local.find(g0);
        if (it0 != global_to_local.end()) {
            int l0 = it0->second;
            result->uvs[g0 * 2 + 0] = island_uvs[l0 * 2 + 0];
            result->uvs[g0 * 2 + 1] = island_uvs[l0 * 2 + 1];
        }
        auto it1 = global_to_local.find(g1);
        if (it1 != global_to_local.end()) {
            int l1 = it1->second;
            result->uvs[g1 * 2 + 0] = island_uvs[l1 * 2 + 0];
            result->uvs[g1 * 2 + 1] = island_uvs[l1 * 2 + 1];
        }
        auto it2 = global_to_local.find(g2);
        if (it2 != global_to_local.end()) {
            int l2 = it2->second;
            result->uvs[g2 * 2 + 0] = island_uvs[l2 * 2 + 0];
            result->uvs[g2 * 2 + 1] = island_uvs[l2 * 2 + 1];
        }
    }
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

    // STEP 1: Build topology
    TopologyInfo* topo = build_topology(mesh);
    if (!topo) {
        fprintf(stderr, "Failed to build topology\n");
        return NULL;
    }
    validate_topology(mesh, topo);

    // STEP 2: Detect seams
    int num_seams = 0;
    int* seam_edges = detect_seams(mesh, topo, params->angle_threshold, &num_seams);

    // STEP 3: Extract islands
    int num_islands = 0;
    int* face_island_ids = extract_islands(mesh, topo, seam_edges, num_seams, &num_islands);

    // STEP 4: Parameterize each island using LSCM
    Mesh* result = allocate_mesh_copy(mesh);
    result->uvs = (float*)calloc(mesh->num_vertices * 2, sizeof(float));
    if (!result->uvs) {
        fprintf(stderr, "Failed to allocate result UVs\n");
        free_topology(topo);
        free(seam_edges);
        free_unwrap_result(NULL);
        return NULL;
    }

    // For temporary container of faces per island
    std::vector<std::vector<int>> island_face_lists(num_islands);
    for (int f = 0; f < mesh->num_triangles; ++f) {
        int id = face_island_ids[f];
        if (id >= 0 && id < num_islands) {
            island_face_lists[id].push_back(f);
        }
    }

    // Buffers and containers used per-island
    for (int island_id = 0; island_id < num_islands; ++island_id) {
        const std::vector<int>& island_faces = island_face_lists[island_id];
        int num_faces_in_island = (int)island_faces.size();
        printf("\nProcessing island %d/%d...\n", island_id + 1, num_islands);
        printf("  %d faces in island\n", num_faces_in_island);

        if (num_faces_in_island == 0) {
            printf("  Skipping (empty)\n");
            continue;
        }
        if (num_faces_in_island < params->min_island_faces) {
            printf("  Skipping (too small)\n");
            continue;
        }

        // Collect unique global vertices used by this island
        std::map<int, int> global_to_local; // global vertex idx -> local idx
        std::vector<int> local_to_global;    // reverse mapping
        local_to_global.reserve(num_faces_in_island * 3);

        for (int fi = 0; fi < num_faces_in_island; ++fi) {
            int fidx = island_faces[fi];
            int tbase = 3 * fidx;
            int g0 = mesh->triangles[tbase + 0];
            int g1 = mesh->triangles[tbase + 1];
            int g2 = mesh->triangles[tbase + 2];
            int globals[3] = {g0, g1, g2};
            for (int k = 0; k < 3; ++k) {
                int gv = globals[k];
                if (global_to_local.find(gv) == global_to_local.end()) {
                    int local_idx = (int)local_to_global.size();
                    global_to_local[gv] = local_idx;
                    local_to_global.push_back(gv);
                }
            }
        }

        int num_local_vertices = (int)local_to_global.size();

        // Build local vertex positions array (float array of size 3 * num_local_vertices)
        float* V_local = (float*)malloc(sizeof(float) * 3 * num_local_vertices);
        if (!V_local) {
            fprintf(stderr, "Failed to allocate V_local\n");
            continue;
        }
        for (int li = 0; li < num_local_vertices; ++li) {
            int gv = local_to_global[li];
            // each mesh vertex is at mesh->vertices[gv*3 + ...]
            V_local[li*3 + 0] = mesh->vertices[gv*3 + 0];
            V_local[li*3 + 1] = mesh->vertices[gv*3 + 1];
            V_local[li*3 + 2] = mesh->vertices[gv*3 + 2];
        }

        // Build local face index array (3 * num_faces_in_island) using local indices
        int* F_local = (int*)malloc(sizeof(int) * 3 * num_faces_in_island);
        if (!F_local) {
            fprintf(stderr, "Failed to allocate F_local\n");
            free(V_local);
            continue;
        }
        for (int fi = 0; fi < num_faces_in_island; ++fi) {
            int fidx = island_faces[fi];
            int tbase = 3 * fidx;
            int g0 = mesh->triangles[tbase + 0];
            int g1 = mesh->triangles[tbase + 1];
            int g2 = mesh->triangles[tbase + 2];

            F_local[3*fi + 0] = global_to_local[g0];
            F_local[3*fi + 1] = global_to_local[g1];
            F_local[3*fi + 2] = global_to_local[g2];
        }

        // Build a temporary Mesh describing the island (vertices pointer + counts).
        // lscm_parameterize expects a Mesh* and an int* face list and the number of faces.
        Mesh local_mesh;
        local_mesh.vertices = V_local;
        local_mesh.num_vertices = num_local_vertices;
        local_mesh.triangles = NULL;          // triangles not used because we pass F_local
        local_mesh.num_triangles = num_faces_in_island;
        local_mesh.uvs = NULL;

        // Call LSCM parameterizer. It returns pointer to floats [u0,v0,u1,v1,...] for local vertices.
        float* UV_local = lscm_parameterize(&local_mesh, F_local, num_faces_in_island);
        if (!UV_local) {
            fprintf(stderr, "LSCM failed for island %d\n", island_id);
            free(V_local);
            free(F_local);
            continue;
        }

        // Copy local UVs into global result->uvs using mapping
        // Prepare an array of face indices for copy function
        std::vector<int> face_indices_vec = island_faces;
        // convert to C array for the helper
        int* face_indices_c = (int*)malloc(sizeof(int) * num_faces_in_island);
        for (int i = 0; i < num_faces_in_island; ++i) face_indices_c[i] = face_indices_vec[i];

        copy_island_uvs(result, UV_local, face_indices_c, num_faces_in_island, global_to_local);

        // Free temporary allocations
        // lscm_parameterize allocated UV_local (contract: caller frees)
        free(UV_local);
        free(V_local);
        free(F_local);
        free(face_indices_c);
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
    if (seam_edges) free(seam_edges);

    printf("\n=== Unwrapping Complete ===\n");
    return result;
}

void free_unwrap_result(UnwrapResult* result) {
    if (!result) return;
    if (result->face_island_ids) free(result->face_island_ids);
    free(result);
}

