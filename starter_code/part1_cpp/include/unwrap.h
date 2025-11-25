/**
 * @file unwrap.h
 * @brief Main UV unwrapping API
 *
 * API SKELETON - YOU IMPLEMENT IN unwrap.cpp, seam_detection.cpp, packing.cpp
 */

#ifndef UNWRAP_H
#define UNWRAP_H

#include "mesh.h"
#include "topology.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Unwrapping parameters
 */
typedef struct {
    float angle_threshold;       /**< Seam detection angle threshold (degrees) */
    int min_island_faces;        /**< Minimum island size (merge smaller islands) */
    int pack_islands;            /**< If true, pack islands into [0,1]² */
    float island_margin;         /**< Spacing between islands (e.g., 0.02) */
} UnwrapParams;

/**
 * @brief Unwrapping result metadata
 */
typedef struct {
    int num_islands;             /**< Number of UV islands */
    int* face_island_ids;        /**< Island ID per face (num_triangles) */
    float avg_stretch;           /**< Average stretch across all triangles */
    float max_stretch;           /**< Maximum stretch */
    float coverage;              /**< Percentage of [0,1]² used */
} UnwrapResult;

/**
 * @brief Main unwrapping function
 *
 * Algorithm:
 * 1. Build mesh topology
 * 2. Detect seams using spanning tree + angular defect
 * 3. Extract UV islands (connected components after seam cuts)
 * 4. Parameterize each island using LSCM
 * 5. Pack islands into [0,1]²
 * 6. Compute quality metrics
 *
 * @param mesh Input mesh
 * @param params Unwrapping parameters
 * @param result_out Output metadata (allocated by function)
 * @return New mesh with UVs, or NULL on error
 * @note Caller must free result mesh and result_out
 *
 * IMPLEMENTATION REQUIRED
 */
Mesh* unwrap_mesh(const Mesh* mesh,
                  const UnwrapParams* params,
                  UnwrapResult** result_out);

/**
 * @brief Detect seams for unwrapping
 *
 * Algorithm:
 * 1. Build dual graph (faces as nodes, shared edges as edges)
 * 2. Compute spanning tree via BFS
 * 3. Mark non-tree edges as seam candidates
 * 4. Refine using angular defect (add edges near high-curvature vertices)
 *
 * @param mesh Input mesh
 * @param topo Topology information
 * @param angle_threshold Angle threshold in degrees
 * @param num_seams_out Output: number of seams detected
 * @return Array of seam edge indices
 * @note Caller must free returned array
 *
 * IMPLEMENTATION REQUIRED - seam_detection.cpp
 */
int* detect_seams(const Mesh* mesh,
                  const TopologyInfo* topo,
                  float angle_threshold,
                  int* num_seams_out);

/**
 * @brief Pack UV islands into [0,1]² texture space
 *
 * Algorithm:
 * 1. Compute bounding box for each island
 * 2. Sort islands by height (descending)
 * 3. Pack using shelf algorithm
 * 4. Scale to fit [0,1]²
 *
 * @param mesh Mesh with UVs (modified in-place)
 * @param result Unwrap result with island IDs
 * @param margin Spacing between islands
 *
 * IMPLEMENTATION REQUIRED - packing.cpp
 */
void pack_uv_islands(Mesh* mesh,
                     const UnwrapResult* result,
                     float margin);

/**
 * @brief Compute quality metrics for UV mapping
 * @param mesh Mesh with UVs
 * @param result Result structure to fill with metrics
 */
void compute_quality_metrics(const Mesh* mesh, UnwrapResult* result);

/**
 * @brief Free unwrap result
 * @param result Result to free
 */
void free_unwrap_result(UnwrapResult* result);

#ifdef __cplusplus
}
#endif

#endif /* UNWRAP_H */
