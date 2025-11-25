/**
 * @file lscm.cpp
 * @brief LSCM (Least Squares Conformal Maps) parameterization
 *
 * SKELETON - YOU IMPLEMENT THIS
 *
 * This is the MOST COMPLEX part of the assignment.
 * Budget 4-5 hours for this file.
 *
 * IMPORTANT:
 * - See reference/lscm_matrix_example.cpp for matrix assembly example
 * - See reference/algorithms.md for mathematical background
 * - Use Eigen library for sparse linear algebra
 *
 * Algorithm:
 * 1. Build local vertex mapping (global → local indices)
 * 2. Assemble LSCM sparse matrix
 * 3. Set boundary conditions (pin 2 vertices)
 * 4. Solve sparse linear system
 * 5. Normalize UVs to [0,1]²
 */

#include "lscm.h"
#include "math_utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <map>
#include <vector>
#include <set>

// Eigen library for sparse matrices
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
// Alternative: #include <Eigen/IterativeLinearSolvers>

int find_boundary_vertices(const Mesh* mesh,
                          const int* face_indices,
                          int num_faces,
                          int** boundary_out) {
    // TODO: Implement
    //
    // Steps:
    // 1. Count how many times each edge appears in the island
    // 2. Edges that appear only once are boundary edges
    // 3. Collect vertices from boundary edges
    //
    // Hint: Use std::map<Edge, int> to count edge usage

    std::set<int> boundary_verts;

    // YOUR CODE HERE

    // Convert to array
    int num_boundary = boundary_verts.size();
    *boundary_out = (int*)malloc(num_boundary * sizeof(int));

    int idx = 0;
    for (int v : boundary_verts) {
        (*boundary_out)[idx++] = v;
    }

    return num_boundary;
}

void normalize_uvs_to_unit_square(float* uvs, int num_verts) {
    if (!uvs || num_verts == 0) return;

    // Find bounding box
    float min_u = FLT_MAX, max_u = -FLT_MAX;
    float min_v = FLT_MAX, max_v = -FLT_MAX;

    for (int i = 0; i < num_verts; i++) {
        float u = uvs[i * 2];
        float v = uvs[i * 2 + 1];

        min_u = min_float(min_u, u);
        max_u = max_float(max_u, u);
        min_v = min_float(min_v, v);
        max_v = max_float(max_v, v);
    }

    float u_range = max_u - min_u;
    float v_range = max_v - min_v;

    if (u_range < 1e-6f) u_range = 1.0f;
    if (v_range < 1e-6f) v_range = 1.0f;

    // Normalize to [0, 1]
    for (int i = 0; i < num_verts; i++) {
        uvs[i * 2] = (uvs[i * 2] - min_u) / u_range;
        uvs[i * 2 + 1] = (uvs[i * 2 + 1] - min_v) / v_range;
    }
}

float* lscm_parameterize(const Mesh* mesh,
                         const int* face_indices,
                         int num_faces) {
    if (!mesh || !face_indices || num_faces == 0) return NULL;

    // TODO: Implement LSCM parameterization
    //
    // This is the CORE algorithm. Take your time.
    //
    // STEP 1: Build local vertex mapping
    //   - Map global vertex indices → local island indices (0, 1, 2, ...)
    //   - Use std::map<int, int> global_to_local
    //   - Use std::vector<int> local_to_global
    //
    // STEP 2: Build LSCM sparse matrix (2n × 2n)
    //   - n = number of vertices in island
    //   - Matrix variables: [u0, v0, u1, v1, ..., u_{n-1}, v_{n-1}]
    //   - For each triangle:
    //       a) Get 3D positions of vertices
    //       b) Project triangle to its plane (create local 2D coords)
    //       c) Compute triangle area (weight)
    //       d) Add LSCM energy terms to matrix
    //   - Use Eigen::Triplet<double> to build matrix
    //   - See reference/lscm_matrix_example.cpp for exact formulas
    //
    // STEP 3: Set boundary conditions
    //   - Find boundary vertices (or pick 2 arbitrary vertices if closed)
    //   - Find 2 vertices far apart
    //   - Pin vertex 1 to (0, 0)
    //   - Pin vertex 2 to (1, 0)
    //   - Modify matrix rows for pinned vertices
    //
    // STEP 4: Solve sparse linear system
    //   - Use Eigen::SparseLU (more robust)
    //   - OR Eigen::ConjugateGradient (faster for large meshes)
    //   - Solve Ax = b
    //
    // STEP 5: Extract and normalize UVs
    //   - Extract u,v from solution vector
    //   - Normalize to [0,1]²
    //
    // PERFORMANCE TARGETS:
    //   - 10,000 vertices with SparseLU: < 5 seconds
    //   - 10,000 vertices with ConjugateGradient: < 2 seconds
    //
    // QUALITY TARGETS:
    //   - Max stretch: < 1.5
    //   - Cylinder test: < 1.2

    printf("LSCM parameterizing %d faces...\n", num_faces);

    // STEP 1: Local vertex mapping
    std::map<int, int> global_to_local;
    std::vector<int> local_to_global;

    // YOUR CODE HERE
    // ...

    int n = local_to_global.size();
    printf("  Island has %d vertices\n", n);

    if (n < 3) {
        fprintf(stderr, "LSCM: Island too small (%d vertices)\n", n);
        return NULL;
    }

    // STEP 2: Build sparse matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> triplets;

    // YOUR CODE HERE
    // For each triangle:
    //   - Get vertices and 3D positions
    //   - Project to triangle plane
    //   - Add LSCM energy terms
    // See reference/lscm_matrix_example.cpp

    // STEP 3: Boundary conditions
    // YOUR CODE HERE
    // Find 2 vertices to pin

    // STEP 4: Solve
    Eigen::SparseMatrix<double> A(2*n, 2*n);
    A.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::VectorXd b = Eigen::VectorXd::Zero(2*n);

    // YOUR CODE HERE
    // Set up solver and solve

    // STEP 5: Extract UVs
    float* uvs = (float*)malloc(n * 2 * sizeof(float));

    // YOUR CODE HERE
    // Extract from solution vector

    normalize_uvs_to_unit_square(uvs, n);

    printf("  LSCM completed\n");
    return uvs;
}
