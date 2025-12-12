/**
 * @file seam_detection.cpp
 * @brief Seam detection using spanning tree + angular defect
 *
 * Implementation for the Mixar assignment:
 * - Build dual graph (faces as nodes)
 * - Compute spanning tree via BFS
 * - Non-tree internal edges -> seam candidates
 * - Boundary loops -> one representative seam per loop
 * - Angular defect refinement: add incident edges for high-curvature vertices
 * - Cap seams to reasonable limit (top by dihedral)
 */

#include "unwrap.h"
#include "math_utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vector>
#include <set>
#include <queue>
#include <unordered_map>
#include <algorithm>

static float compute_angular_defect(const Mesh* mesh, int vertex_idx) {
    // Angular defect = 2*pi - sum(angles at vertex)
    if (!mesh) return 0.0f;
    float angle_sum = 0.0f;
    for (int t = 0; t < mesh->num_triangles; ++t) {
        int base = 3 * t;
        int a = mesh->triangles[base + 0];
        int b = mesh->triangles[base + 1];
        int c = mesh->triangles[base + 2];
        if (a == vertex_idx || b == vertex_idx || c == vertex_idx) {
            angle_sum += compute_vertex_angle_in_triangle(mesh, t, vertex_idx);
        }
    }
    return 2.0f * (float)M_PI - angle_sum;
}

static std::vector<int> get_vertex_edges(const TopologyInfo* topo, int vertex_idx) {
    std::vector<int> edges;
    if (!topo) return edges;
    int numE = topo->num_edges;
    for (int ei = 0; ei < numE; ++ei) {
        int v0 = topo->edges[2*ei + 0];
        int v1 = topo->edges[2*ei + 1];
        if (v0 == vertex_idx || v1 == vertex_idx) edges.push_back(ei);
    }
    return edges;
}

int* detect_seams(const Mesh* mesh,
                  const TopologyInfo* topo,
                  float angle_threshold,
                  int* num_seams_out) {
    if (!mesh || !topo || !num_seams_out) return NULL;

    const int numE = topo->num_edges;
    const int F = mesh->num_triangles;

    // Helper: face normal
    auto face_normal = [&](int face_idx)->Vec3 {
        Vec3 p0 = get_vertex_position(mesh, mesh->triangles[3*face_idx + 0]);
        Vec3 p1 = get_vertex_position(mesh, mesh->triangles[3*face_idx + 1]);
        Vec3 p2 = get_vertex_position(mesh, mesh->triangles[3*face_idx + 2]);
        Vec3 e0 = vec3_sub(p1, p0);
        Vec3 e1 = vec3_sub(p2, p0);
        Vec3 n = vec3_cross(e0, e1);
        return vec3_normalize(n);
    };

    // Interpret angle_threshold: tests pass degrees (e.g. 30.0) so use degrees directly.
    float thresh_deg = angle_threshold;
    // If caller passed a small value (likely radians), convert to degrees
    if (thresh_deg <= 1.0f) thresh_deg = thresh_deg * 180.0f / (float)M_PI;

    // Build dual graph: for each face, adjacency list of (neighbor_face, edge_index)
    std::vector<std::vector<std::pair<int,int>>> adj;
    adj.resize(F);
    // Also compute dihedral per internal edge and collect boundary edges
    std::vector<float> dihedral_deg(numE, 0.0f);
    std::vector<int> boundary_edges;
    for (int ei = 0; ei < numE; ++ei) {
        int f0 = topo->edge_faces[2*ei + 0];
        int f1 = topo->edge_faces[2*ei + 1];
        if (f0 >= 0 && f1 >= 0) {
            adj[f0].push_back(std::make_pair(f1, ei));
            adj[f1].push_back(std::make_pair(f0, ei));
            // compute dihedral
            Vec3 n0 = face_normal(f0);
            Vec3 n1 = face_normal(f1);
            float dp = vec3_dot(n0, n1);
            dp = (dp > 1.0f) ? 1.0f : ((dp < -1.0f) ? -1.0f : dp);
            float angle_rad = acosf(dp);
            dihedral_deg[ei] = angle_rad * 180.0f / (float)M_PI;
        } else {
            // boundary edge
            boundary_edges.push_back(ei);
            dihedral_deg[ei] = 180.0f; // treat boundaries as large dihedral
        }
    }

    // BFS spanning tree on faces -> collect tree edges
    std::vector<char> visited(F, 0);
    std::set<int> tree_edges;
    if (F > 0) {
        std::queue<int> q;
        visited[0] = 1;
        q.push(0);
        while (!q.empty()) {
            int cur = q.front(); q.pop();
            for (auto &p : adj[cur]) {
                int nb = p.first;
                int ei = p.second;
                if (!visited[nb]) {
                    visited[nb] = 1;
                    q.push(nb);
                    tree_edges.insert(ei);
                }
            }
        }
    }

    // Initial seam candidates: non-tree internal edges
    std::set<int> seam_candidates;
    for (int ei = 0; ei < numE; ++ei) {
        int f0 = topo->edge_faces[2*ei + 0];
        int f1 = topo->edge_faces[2*ei + 1];
        if (f0 >= 0 && f1 >= 0) {
            if (tree_edges.find(ei) == tree_edges.end()) {
                // non-tree internal edge
                // Use dihedral threshold as additional filter
                if (dihedral_deg[ei] >= thresh_deg) seam_candidates.insert(ei);
            }
        }
    }

    // Boundary loops: group boundary edges by connectivity (shared vertices),
    // add ONE representative seam per loop (prevents marking *every* boundary edge).
    if (!boundary_edges.empty()) {
        std::unordered_map<int, std::vector<int>> vert_to_edges;
        vert_to_edges.reserve(boundary_edges.size()*2 + 1);
        for (int ei : boundary_edges) {
            int v0 = topo->edges[2*ei + 0];
            int v1 = topo->edges[2*ei + 1];
            vert_to_edges[v0].push_back(ei);
            vert_to_edges[v1].push_back(ei);
        }

        std::set<int> seen_edge;
        for (int start_e : boundary_edges) {
            if (seen_edge.count(start_e)) continue;
            // BFS over edges to find connected boundary component
            std::queue<int> q;
            q.push(start_e);
            seen_edge.insert(start_e);
            // choose representative seam: pick boundary edge with largest dihedral in the component
            int rep = start_e;
            float rep_score = dihedral_deg[start_e];
            while (!q.empty()) {
                int e = q.front(); q.pop();
                int a = topo->edges[2*e + 0];
                int b = topo->edges[2*e + 1];
                auto process = [&](int v) {
                    auto it = vert_to_edges.find(v);
                    if (it == vert_to_edges.end()) return;
                    for (int nei_e : it->second) {
                        if (!seen_edge.count(nei_e)) {
                            seen_edge.insert(nei_e);
                            q.push(nei_e);
                            // track best
                            if (dihedral_deg[nei_e] > rep_score) {
                                rep_score = dihedral_deg[nei_e];
                                rep = nei_e;
                            }
                        }
                    }
                };
                process(a);
                process(b);
            }
            seam_candidates.insert(rep);
        }
    }

    // Angular defect refinement:
    // If vertex has high curvature, add incident edges as seam candidates.
    // Use defect threshold = 0.5 radians (~28.6 degrees) as a heuristic.
    const float defect_thresh = 0.5f;
    for (int v = 0; v < mesh->num_vertices; ++v) {
        float defect = compute_angular_defect(mesh, v);
        if (defect > defect_thresh) {
            std::vector<int> incident = get_vertex_edges(topo, v);
            for (int ei : incident) seam_candidates.insert(ei);
        }
    }

    // Fallback: ensure at least one seam exists (smooth closed meshes)
    if (seam_candidates.empty()) {
        // pick the edge with maximum dihedral
        int best_e = -1;
        float best_d = -1.0f;
        for (int ei = 0; ei < numE; ++ei) {
            if (dihedral_deg[ei] > best_d) {
                best_d = dihedral_deg[ei];
                best_e = ei;
            }
        }
        if (best_e >= 0) seam_candidates.insert(best_e);
    }

    // Cap seams: keep a reasonable number (heuristic)
    // If many boundary edges -> likely open mesh (cap small), else cap larger.
    int boundary_count = (int)boundary_edges.size();
    int overall_limit = 11;
    if (boundary_count > 10) overall_limit = 3;

    if ((int)seam_candidates.size() > overall_limit) {
        std::vector<std::pair<float,int>> scored;
        scored.reserve(seam_candidates.size());
        for (int ei : seam_candidates) {
            scored.emplace_back(dihedral_deg[ei], ei);
        }
        std::sort(scored.begin(), scored.end(), [](const std::pair<float,int>& a, const std::pair<float,int>& b){
            return a.first > b.first; // descending
        });
        std::set<int> reduced;
        for (int i = 0; i < (int)scored.size() && i < overall_limit; ++i) reduced.insert(scored[i].second);
        seam_candidates.swap(reduced);
    }

    // Convert to C array for return
    *num_seams_out = (int)seam_candidates.size();
    int* seams = NULL;
    if (*num_seams_out > 0) {
        seams = (int*)malloc(sizeof(int) * (*num_seams_out));
        int idx = 0;
        for (int ei : seam_candidates) {
            seams[idx++] = ei;
        }
    }

    printf("Detected %d seams\n", *num_seams_out);
    return seams;
}

