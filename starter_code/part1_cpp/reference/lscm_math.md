# LSCM (Least Squares Conformal Maps) - Mathematical Background

## Overview

LSCM is a powerful method for computing conformal parameterization of 3D surfaces. This document provides the mathematical foundation you need to implement the algorithm.

## 1. Conformal Maps

A mapping `φ: S → R²` from a surface S to the plane is **conformal** if it preserves angles locally. Mathematically, this means the Jacobian of φ at each point is a similarity transformation (rotation + uniform scaling).

### Why Conformal?

- Preserves shape locally (no shearing)
- Minimizes distortion for texture mapping
- Optimal for UV unwrapping in graphics

## 2. The LSCM Energy Functional

For a triangle mesh, LSCM minimizes the conformal energy:

```
E(u,v) = ∑ₜ Area(t) · ||∇u - R₉₀°(∇v)||²
```

Where:
- `t` ranges over all triangles
- `u, v` are the UV coordinates for each vertex
- `∇u, ∇v` are gradients of u and v over triangle t
- `R₉₀°` is the 90° counterclockwise rotation operator
- For a vector `∇v = [vₓ, vᵧ]`, the rotation is `R₉₀°(∇v) = [-vᵧ, vₓ]`

### Physical Interpretation

The energy measures how much the mapping deviates from being conformal. When E = 0, the mapping is perfectly conformal (angles preserved).

## 3. Discretization for Triangle Meshes

For a triangle with vertices `(x₀,y₀,z₀)`, `(x₁,y₁,z₁)`, `(x₂,y₂,z₂)` and UV coordinates `(u₀,v₀)`, `(u₁,v₁)`, `(u₂,v₂)`:

### Step 1: Compute Local 2D Basis

Project the triangle to a local 2D coordinate system:

```
X₁ = (x₁ - x₀, y₁ - y₀, z₁ - z₀)
X₂ = (x₂ - x₀, y₂ - y₀, z₂ - z₀)

e₁ = X₁ / ||X₁||  (unit vector along first edge)
e₂ = (X₂ - (X₂·e₁)e₁) / ||X₂ - (X₂·e₁)e₁||  (orthogonal unit vector)
```

Local 2D coordinates:
```
x'₀ = 0,  y'₀ = 0
x'₁ = ||X₁||,  y'₁ = 0
x'₂ = X₂·e₁,  y'₂ = X₂·e₂
```

### Step 2: Compute Gradients

For a triangle in 2D, the gradient of any function is constant:

```
∇u = (1/2A) [(y'₁ - y'₂)(u₀) + (y'₂ - y'₀)(u₁) + (y'₀ - y'₁)(u₂),
              (x'₂ - x'₁)(u₀) + (x'₀ - x'₂)(u₁) + (x'₁ - x'₀)(u₂)]

∇v = (1/2A) [(y'₁ - y'₂)(v₀) + (y'₂ - y'₀)(v₁) + (y'₀ - y'₁)(v₂),
              (x'₂ - x'₁)(v₀) + (x'₀ - x'₂)(v₁) + (x'₁ - x'₀)(v₂)]
```

Where `A` is the triangle area: `A = 0.5 * |x'₁ * y'₂ - x'₂ * y'₁|`

### Step 3: Conformal Energy for One Triangle

```
Eₜ = A · [((∇u)ₓ - (∇v)ᵧ)² + ((∇u)ᵧ + (∇v)ₓ)²]
```

This directly relates to the Cauchy-Riemann equations. For a conformal map:
- `∂u/∂x = ∂v/∂y` (satisfied when `(∇u)ₓ - (∇v)ᵧ = 0`)
- `∂u/∂y = -∂v/∂x` (satisfied when `(∇u)ᵧ + (∇v)ₓ = 0`)

The energy is zero when both conditions are met. This simplifies to a quadratic form in the UV coordinates.

## 4. Matrix Assembly

The total energy is:

```
E = ∑ₜ Eₜ = UᵀMU
```

Where:
- `U = [u₀, v₀, u₁, v₁, ..., uₙ, vₙ]ᵀ` is the vector of all UV coordinates
- `M` is a sparse (2n × 2n) matrix assembled from triangle contributions

### Matrix Construction

For each triangle `t` with vertices `i, j, k`:

1. Compute local 2D coordinates `(x'ᵢ, y'ᵢ)` for each vertex
2. Compute area `A`
3. Compute coefficient matrix for this triangle:

```
For vertex pairs (p,q) in {i,j,k}:
  Compute contributions to M[2p:2p+2, 2q:2q+2]
  using the conformal energy formula
```

The matrix M is symmetric and sparse (only entries corresponding to edge-connected vertices are non-zero).

## 5. Boundary Conditions

The LSCM energy has a null space (translation, rotation, uniform scaling). To get a unique solution, we must fix some DOFs:

### Standard Approach: Pin 2 Vertices

1. Find two boundary vertices `v₁, v₂` (or any two vertices for closed surfaces)
2. Fix their UV coordinates:
   ```
   u₁ = 0,  v₁ = 0
   u₂ = 1,  v₂ = 0
   ```

3. Solve for remaining unknowns

### Implementation

Partition variables into free and fixed:
```
U = [U_free; U_fixed]
```

The energy becomes:
```
E = U_freeᵀ M₁₁ U_free + 2U_freeᵀ M₁₂ U_fixed + const
```

Minimizing ∂E/∂U_free = 0 gives:
```
M₁₁ U_free = -M₁₂ U_fixed
```

## 6. Solving the Linear System

The system `M₁₁ U_free = b` is sparse, symmetric, and positive semi-definite.

### Using Eigen Library

```cpp
#include <Eigen/Sparse>

Eigen::SparseMatrix<double> M(n, n);
Eigen::VectorXd b(n);

// ... assemble M and b ...

// Solve using conjugate gradient
Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
solver.compute(M);
Eigen::VectorXd U = solver.solve(b);
```

### Alternative Solvers

- `SimplicialLLT` (faster for small systems)
- `SparseLU` (more robust, slower)
- `ConjugateGradient` (best for large systems)

## 7. Post-processing

After solving, normalize UVs to [0,1]²:

```cpp
// Find bounding box
float u_min = min(all u values);
float u_max = max(all u values);
float v_min = min(all v values);
float v_max = max(all v values);

// Normalize
for each vertex i:
    u[i] = (u[i] - u_min) / (u_max - u_min);
    v[i] = (v[i] - v_min) / (v_max - v_min);
```

## 8. Implementation Pseudocode

```
function lscm_unwrap(mesh, boundary_faces):
    // 1. Build local vertex mapping
    vertices = unique_vertices(boundary_faces)
    n = vertices.size()
    global_to_local = build_mapping(vertices)

    // 2. Find boundary vertices to pin
    boundary_verts = find_boundary(mesh, boundary_faces)
    pin1 = boundary_verts[0]
    pin2 = boundary_verts[boundary_verts.size()/2]

    // 3. Assemble LSCM matrix
    M = SparseMatrix(2n, 2n)
    for each face in boundary_faces:
        v0, v1, v2 = face.vertices
        add_triangle_to_matrix(M, v0, v1, v2)

    // 4. Apply boundary conditions
    M_reduced, b = fix_vertices(M, pin1, pin2)

    // 5. Solve
    solver = ConjugateGradient()
    solver.compute(M_reduced)
    UV_free = solver.solve(b)

    // 6. Reconstruct full solution
    UV = reconstruct(UV_free, pin1, pin2)

    // 7. Normalize
    normalize_uvs(UV)

    return UV
```

## 9. Common Pitfalls

### Singular Matrix
- Forgot to apply boundary conditions (pin vertices)
- Disconnected components in mesh
- Solution: Always pin at least 2 vertices with distinct positions

### Poor Quality UVs
- Boundary conditions too close together
- High distortion at seams
- Solution: Choose well-separated boundary vertices

### Numerical Instability
- Very small triangles (near-zero area)
- Degenerate triangles (collinear vertices)
- Solution: Add regularization term: `M + εI`

## 10. Advanced Topics

### Free Boundary LSCM

Instead of pinning vertices, constrain only translation/rotation:
- Minimize `||UV - R*UV₀ - t||²` where R is rotation
- Requires iterative solver (e.g., alternating optimization)

### Multi-resolution LSCM

1. Simplify mesh
2. Solve LSCM on coarse mesh
3. Upsample UVs
4. Refine with fine-level LSCM

### Seamless Parameterization

For closed surfaces without cutting:
- Use multiple charts
- Optimize transition functions between charts
- Advanced topic beyond this assignment

## References

1. Lévy, B., et al. (2002). "Least Squares Conformal Maps for Automatic Texture Atlas Generation." ACM SIGGRAPH 2002.

2. Sheffer, A., et al. (2006). "Mesh Parameterization Methods and Their Applications." Foundations and Trends in Computer Graphics and Vision.

3. Eigen library documentation: https://eigen.tuxfamily.org/

## Quick Reference

**Key equations:**
```
Energy:     E = ∑ₜ Area(t) · ||(∇u)ᵀ - (∇v)ᴶ||²
Matrix:     M₁₁ U_free = -M₁₂ U_fixed
Normalize:  u' = (u - u_min)/(u_max - u_min)
```

**Algorithm steps:**
1. Build local vertex mapping (global → local indices)
2. Find 2 boundary vertices to pin
3. Assemble sparse LSCM matrix M
4. Apply boundary conditions (reduce system)
5. Solve M₁₁ U_free = b using CG solver
6. Reconstruct full UV vector
7. Normalize to [0,1]²

Good luck with your implementation!
