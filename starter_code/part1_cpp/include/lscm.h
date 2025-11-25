/**
 * @file lscm.h
 * @brief LSCM (Least Squares Conformal Maps) parameterization
 *
 * API SKELETON - YOU IMPLEMENT IN lscm.cpp
 */

#ifndef LSCM_H
#define LSCM_H

#include "mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Parameterize a UV island using LSCM
 *
 * Algorithm:
 * 1. Build local vertex mapping (global indices → local indices)
 * 2. Assemble LSCM sparse matrix (2n × 2n)
 *    - For each triangle, add energy contribution
 *    - Energy: ||∇u - R_90°(∇v)||²
 * 3. Set boundary conditions (pin 2 vertices to prevent degeneracy)
 * 4. Solve sparse linear system
 * 5. Normalize UVs to [0,1]²
 *
 * @param mesh Input mesh
 * @param face_indices Indices of faces in this island
 * @param num_faces Number of faces in island
 * @return Array of UVs [u,v, u,v, ...] for vertices in island
 * @note Caller must free returned array
 *
 * IMPLEMENTATION REQUIRED
 *
 * IMPORTANT: See reference/lscm_matrix_example.cpp for matrix assembly example
 *
 * Dependencies:
 * - Eigen library for sparse matrices (include in third_party/ or update CMakeLists.txt if needed)
 * - Use Eigen::SparseLU or Eigen::ConjugateGradient for solving the linear system
 */
float* lscm_parameterize(const Mesh* mesh,
                         const int* face_indices,
                         int num_faces);

/**
 * @brief Helper: Find boundary vertices in an island
 * @param mesh Input mesh
 * @param face_indices Faces in island
 * @param num_faces Number of faces
 * @param boundary_out Output array of boundary vertex indices
 * @return Number of boundary vertices
 */
int find_boundary_vertices(const Mesh* mesh,
                          const int* face_indices,
                          int num_faces,
                          int** boundary_out);

/**
 * @brief Helper: Normalize UVs to unit square [0,1]²
 * @param uvs UV array to normalize (modified in-place)
 * @param num_verts Number of vertices
 */
void normalize_uvs_to_unit_square(float* uvs, int num_verts);

#ifdef __cplusplus
}
#endif

#endif /* LSCM_H */
