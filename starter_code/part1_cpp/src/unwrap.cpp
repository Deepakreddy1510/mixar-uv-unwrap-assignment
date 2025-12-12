/**
 * @file unwrap.cpp
 * @brief Main UV unwrapping orchestrator
 *
 * Implementation filling missing functions:
 * - extract_islands
 * - copy_island_uvs
 * - island LSCM call / local mesh assembly
 */

#include "unwrap.h"
#include "lscm.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <map>
#include <queue>
#include <string.h>

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
    int F = mesh->num_triangles;
    int E = topo->num_edges;

    // Build seam set for quick lookup
    std::set<int> seam_set;
    for (int i = 0; i < num_seams; ++i) {
        seam_set.insert(seam_edges[i]);
    }

    // Build face adjacency via non-seam internal edges
    std::vector<std::vector<int>> adj(F);
    for (int ei = 0; ei < E; ++ei) {
        // Skip seam edges entirely
        if (seam_set.count(ei)) continue;

        int f0 = topo->edge_faces[2*ei + 0];
        int f1 = topo->edge_faces[2*ei + 1];

        // If both faces valid and not boundary, connect them
        if (f0 >= 0 && f1 >= 0) {
            adj[f0].push_back(f1);
            adj[f1].push_back(f0);
        }
    }

    // Prepare island ids, -1 = unvisited
    int* face_island_ids = (int*)malloc(F * sizeof(int));
    if (!face_island_ids) return NULL;
    for (int i = 0; i < F; ++i) face_island_ids[i] = -1;

    int island_id = 0;
    std::vector<int> queue_buf;
    queue_buf.reserve(F);

    for (int f = 0; f < F; ++f) {
        if (face_island_ids[f] != -1) continue; // already assigned

        // BFS from face f
        std::queue<int> q;
        q.push(f);
        face_island_ids[f] = island_id;

        while (!q.empty()) {
            int cur = q.front(); q.pop();
            for (int nei : adj[cur]) {
                if (face_island_ids[nei] == -1) {
                    face_island_ids[nei] = island_id;
                    q.push(nei);
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
 * island_uvs: array of size (2 * num_local_vertices) (float)
 * face_indices: pointer to faces included in island (face indices in original mesh)
 * global_to_local: map from global vertex index -> local vertex index
 */
static void copy_island_uvs(Mesh* result,
                           const float* island_uvs,
                           const int* face_indices,
                           int num_faces,
                           const std::map<int, int>& global_to_local) {
    if (!result || !island_uvs || !face_indices) return;

    // For each global vertex in map, copy its local UV into global result
    // But in some meshes vertices are shared across faces; the mapping covers the island vertices only.
    for (const auto& kv : global_to_local) {
        int global_idx = kv.first;
        int local_idx = kv.second;
        // Each UV is two floats
        result->uvs[2 * global_idx + 0] = island_uvs[2 * local_idx + 0];
        result->uvs[2 * global_idx + 1] = island_uvs[2 * local_idx + 1];
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
    if (!seam_edges && num_seams != 0) {
        // detect_seams should return array when num_seams > 0
        fprintf(stderr, "detect_seams returned null but num_seams=%d\n", num_seams);
        // continue with zero seams
        num_seams = 0;
    }

    // STEP 3: Extract islands
    int num_islands = 0;
    int* face_island_ids = extract_islands(mesh, topo, seam_edges ? seam_edges : (int*)nullptr, num_seams, &num_islands);
    if (!face_island_ids) {
        fprintf(stderr, "Failed to extract islands\n");
        free_topology(topo);
        if (seam_edges) free(seam_edges);
        return NULL;
    }

    // STEP 4: Parameterize each island using LSCM
    Mesh* result = allocate_mesh_copy(mesh);
    // allocate UVs (2 floats per vertex)
    result->uvs = (float*)calloc(mesh->num_vertices * 2, sizeof(float));
    if (!result->uvs) {
        fprintf(stderr, "Failed to allocate result UVs\n");
        free(face_island_ids);
        free_topology(topo);
        if (seam_edges) free(seam_edges);
        free(result);
        return NULL;
    }

    // For each island: collect faces, vertices, build local mesh, call lscm_parameterize, copy back
    for (int island_id = 0; island_id < num_islands; ++island_id) {
        printf("\nProcessing island %d/%d...\n", island_id + 1, num_islands);

        // collect faces
        std::vector<int> island_faces;
        for (int f = 0; f < mesh->num_triangles; ++f) {
            if (face_island_ids[f] == island_id) island_faces.push_back(f);
        }

        int num_faces_in_island = (int)island_faces.size();
        printf("  %d faces in island\n", num_faces_in_island);

        if (num_faces_in_island < params->min_island_faces) {
            printf("  Skipping (too small)\n");
            continue;
        }

        // Build set of unique vertices in this island and mapping global->local
        std::map<int,int> global_to_local;
        std::vector<int> local_to_global;
        local_to_global.reserve(num_faces_in_island * 3);

        for (int fi = 0; fi < num_faces_in_island; ++fi) {
            int fidx = island_faces[fi];
            int t0 = mesh->triangles[3*fidx + 0];
            int t1 = mesh->triangles[3*fidx + 1];
            int t2 = mesh->triangles[3*fidx + 2];
            int tri_v[3] = {t0,t1,t2};
            for (int k = 0; k < 3; ++k) {
                int gv = tri_v[k];
                if (global_to_local.find(gv) == global_to_local.end()) {
                    int local_idx = (int)local_to_global.size();
                    global_to_local[gv] = local_idx;
                    local_to_global.push_back(gv);
                }
            }
        }

        int V_local = (int)local_to_global.size();
        int Nt = num_faces_in_island;

        if (V_local == 0 || Nt == 0) {
            printf("  Island has 0 vertices or 0 faces -> skipping\n");
            continue;
        }

        // Build local Mesh structure for LSCM
        Mesh local_mesh;
        // allocate local vertex positions (3 * V_local floats)
        float* local_vertices = (float*)malloc(sizeof(float) * 3 * V_local);
        if (!local_vertices) {
            fprintf(stderr, "  Failed to allocate local vertices\n");
            continue;
        }
        for (int i = 0; i < V_local; ++i) {
            int gv = local_to_global[i];
            local_vertices[3*i + 0] = mesh->vertices[3*gv + 0];
            local_vertices[3*i + 1] = mesh->vertices[3*gv + 1];
            local_vertices[3*i + 2] = mesh->vertices[3*gv + 2];
        }

        // allocate local triangles (3 * Nt ints)
        int* local_tris = (int*)malloc(sizeof(int) * 3 * Nt);
        if (!local_tris) {
            fprintf(stderr, "  Failed to allocate local triangles\n");
            free(local_vertices);
            continue;
        }
        for (int fi = 0; fi < Nt; ++fi) {
            int fidx = island_faces[fi];
            int gv0 = mesh->triangles[3*fidx + 0];
            int gv1 = mesh->triangles[3*fidx + 1];
            int gv2 = mesh->triangles[3*fidx + 2];
            local_tris[3*fi + 0] = global_to_local[gv0];
            local_tris[3*fi + 1] = global_to_local[gv1];
            local_tris[3*fi + 2] = global_to_local[gv2];
        }

        // Fill local_mesh struct
        local_mesh.vertices = local_vertices;
        local_mesh.num_vertices = V_local;
        local_mesh.triangles = local_tris;
        local_mesh.num_triangles = Nt;
        local_mesh.uvs = NULL; // not needed for input

        // Call LSCM: signature: float* lscm_parameterize(const Mesh* mesh, const int* F_local, int Nt)
        printf("LSCM parameterizing %d faces...\n", Nt);
        float* island_uvs = lscm_parameterize(&local_mesh, local_tris, Nt);
        if (!island_uvs) {
            printf("  LSCM failed for island %d\n", island_id);
            free(local_vertices);
            free(local_tris);
            continue;
        }

        // Copy island UVs into global result using mapping
        // island_uvs is 2 * V_local floats: uv for local vertex i => island_uvs[2*i + 0..1]
        copy_island_uvs(result, island_uvs, island_faces.empty() ? nullptr : &island_faces[0], num_faces_in_island, global_to_local);

        // free island_uvs (LSCM returns allocated float* per its contract)
        free(island_uvs);

        // cleanup local buffers
        free(local_vertices);
        free(local_tris);
    } // end for islands

    // STEP 5: Pack islands if requested
    if (params->pack_islands) {
        UnwrapResult temp_result;
        temp_result.num_islands = num_islands;
        temp_result.face_island_ids = face_island_ids;
        pack_uv_islands(result, &temp_result, params->island_margin);
    }

    // STEP 6: Compute quality metrics
    UnwrapResult* result_data = (UnwrapResult*)malloc(sizeof(UnwrapResult));
    if (!result_data) {
        fprintf(stderr, "Failed to allocate result metadata\n");
        // still return result, but null out result_out
        *result_out = NULL;
    } else {
        result_data->num_islands = num_islands;
        result_data->face_island_ids = face_island_ids;
        compute_quality_metrics(result, result_data);
        *result_out = result_data;
    }

    // Cleanup topology and seam array (face_island_ids is now owned by result_data)
    free_topology(topo);
    if (seam_edges) free(seam_edges);

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

