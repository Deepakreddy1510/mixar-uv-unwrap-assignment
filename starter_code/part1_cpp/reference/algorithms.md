# UV Unwrapping Algorithms Reference

This document provides detailed algorithm descriptions for implementing the UV unwrapping system.

## Table of Contents
1. [Topology Building](#topology-building)
2. [Seam Detection](#seam-detection)
3. [LSCM Parameterization](#lscm-parameterization)
4. [Island Packing](#island-packing)

---

## Topology Building

### Goal
Build a data structure that stores:
- All unique edges in the mesh
- For each edge, the 1 or 2 adjacent faces

### Algorithm

**Input:** Triangle mesh with V vertices and F faces
**Output:** Edge list with adjacency information

```
1. Initialize edge_map as std::map<Edge, EdgeInfo>

2. For each triangle (face) in mesh:
     Get vertices v0, v1, v2

     Add edge (v0, v1) to edge_map with face index
     Add edge (v1, v2) to edge_map with face index
     Add edge (v2, v0) to edge_map with face index

3. Convert edge_map to arrays:
     edges[] = array of (v0, v1) pairs
     edge_faces[] = array of (face0, face1) pairs
                    face1 = -1 for boundary edges

4. Validate using Euler characteristic:
     V - E + F = 2 (for closed meshes)
```

### Edge Uniqueness

Always store edges with smaller vertex index first:
```
Edge(a, b):
    if a < b: return (a, b)
    else: return (b, a)
```

This ensures each edge appears only once in the map.

### Expected Results

| Mesh | Vertices | Edges | Faces | Euler (V-E+F) |
|------|----------|-------|-------|---------------|
| Cube | 8 | 18 | 12 | 2 |
| Sphere | 42 | 120 | 80 | 2 |
| Cylinder (open) | 50 | 96 | 48 | 2 |

---

## Seam Detection

### Goal
Find edges to cut so the mesh can be flattened into 2D without overlaps.

### Algorithm Overview

Uses **spanning tree** on the dual graph plus **angular defect** refinement.

### Dual Graph

- **Nodes:** Faces of the mesh
- **Edges:** Shared edges between faces

```
For each interior edge (edge with 2 adjacent faces):
    face0 and face1 are connected in dual graph
```

### Spanning Tree

A spanning tree visits all faces using F-1 edges.

**Algorithm (BFS):**
```
1. visited[] = false for all faces
2. queue = [face_0]
3. visited[0] = true
4. tree_edges = empty set

5. While queue not empty:
     current_face = queue.pop()

     For each edge adjacent to current_face:
         neighbor_face = other face sharing this edge

         If not visited[neighbor_face]:
             visited[neighbor_face] = true
             tree_edges.add(edge)
             queue.push(neighbor_face)

6. Return tree_edges
```

### Initial Seam Candidates

All edges NOT in the spanning tree are initial seam candidates.

**Reasoning:** Cutting non-tree edges makes the dual graph a tree, which means the mesh becomes topologically a disk (can be flattened).

### Angular Defect Refinement

Some vertices have high curvature and should be on seams.

**Angular defect at vertex v:**
```
defect = 2π - Σ(angles at v in all incident triangles)

For flat surface: defect ≈ 0
For corner (cube): defect > 0
For saddle: defect < 0
```

**Refinement:**
```
For each vertex v:
    defect = compute_angular_defect(v)

    If |defect| > 0.5 radians:  # High curvature
        For each edge incident to v:
            If edge not in spanning tree:
                Add edge to seam candidates
```

### Expected Seam Counts

| Mesh | Expected Seams | Acceptable Range |
|------|----------------|------------------|
| Cube | 7-9 | ±2 |
| Sphere | 1-3 | ±2 |
| Cylinder | 1-2 | ±1 |
| Torus | 2-3 | ±1 |

---

## LSCM Parameterization

### Goal
Find UV coordinates that minimize angle distortion (conformal mapping).

### Mathematical Background

**Conformal mapping:** Preserves angles between curves.

**LSCM Energy:**
```
E = Σ (over triangles) ||∇u - R_90°(∇v)||²

Where:
- ∇u, ∇v are gradients of U and V in 3D space
- R_90° is 90-degree counterclockwise rotation

In component form (using Cauchy-Riemann equations):
E = Σ Area(t) · [((∂u/∂x) - (∂v/∂y))² + ((∂u/∂y) + (∂v/∂x))²]
```

### Discrete Formulation

For each triangle, we build a contribution to a sparse linear system.

**Per-Triangle Contribution:**

Given triangle with vertices p0, p1, p2 in 3D:

1. **Project to triangle's plane:**
   ```
   e1 = p1 - p0
   e2 = p2 - p0
   normal = normalize(cross(e1, e2))

   u_axis = normalize(e1)
   v_axis = cross(normal, u_axis)

   Local 2D coordinates:
   q0 = (0, 0)
   q1 = (dot(e1, u_axis), dot(e1, v_axis))
   q2 = (dot(e2, u_axis), dot(e2, v_axis))
   ```

2. **Triangle area (weight):**
   ```
   area = 0.5 * |q1.x * q2.y - q1.y * q2.x|
   ```

3. **Add energy terms to matrix:**
   See `lscm_matrix_example.cpp` for exact formulas.

   For edge v0→v1:
   ```
   dx = q1.x - q0.x
   dy = q1.y - q0.y

   Matrix[2*v0,   2*v1]   += area * dx
   Matrix[2*v0,   2*v1+1] += area * dy
   Matrix[2*v0+1, 2*v1]   += area * dy
   Matrix[2*v0+1, 2*v1+1] += area * (-dx)

   Matrix[2*v0,   2*v0]   += -area * dx
   Matrix[2*v0,   2*v0+1] += -area * dy
   Matrix[2*v0+1, 2*v0]   += -area * dy
   Matrix[2*v0+1, 2*v0+1] += -area * (-dx)
   ```

   Similarly for edges v1→v2 and v2→v0.

### Boundary Conditions

The system is under-constrained. We must pin 2 vertices to fix translation and rotation.

**Strategy:**
1. Find boundary vertices (edges that appear only once)
2. Find 2 boundary vertices far apart
3. Pin vertex 1 to (0, 0)
4. Pin vertex 2 to (1, 0)

**Implementation:**
```
For pinned vertex at index i with target (u, v):
    Clear row 2*i and row 2*i+1
    Set Matrix[2*i, 2*i] = 1.0
    Set Matrix[2*i+1, 2*i+1] = 1.0
    Set b[2*i] = u
    Set b[2*i+1] = v
```

### Solving

Use Eigen library:

**Option 1: Direct solver (more robust)**
```cpp
Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
solver.compute(A);
Eigen::VectorXd x = solver.solve(b);
```

**Option 2: Iterative solver (faster for large meshes)**
```cpp
Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
solver.setMaxIterations(1000);
solver.setTolerance(1e-6);
solver.compute(A);
Eigen::VectorXd x = solver.solve(b);
```

### Normalization

After solving, normalize UVs to [0,1]²:
```
1. Find min/max of U and V
2. For each UV:
     u_normalized = (u - min_u) / (max_u - min_u)
     v_normalized = (v - min_v) / (max_v - min_v)
```

### Performance Targets

| Mesh Size | SparseLU | ConjugateGradient |
|-----------|----------|-------------------|
| 1k vertices | < 0.5s | < 0.2s |
| 10k vertices | < 5s | < 2s |

### Quality Targets

- Max stretch: < 1.5
- Cylinder test: < 1.2

---

## Island Packing

### Goal
Arrange UV islands efficiently in [0,1]² texture space.

### Shelf Packing Algorithm

Simple and effective algorithm achieving > 60% coverage.

**Algorithm:**
```
1. Compute bounding box for each island:
     For each island:
         min_u, max_u, min_v, max_v
         width = max_u - min_u
         height = max_v - min_v

2. Sort islands by height (descending):
     Larger islands first

3. Pack using shelves:
     shelves = [(x=0, y=0, height=first_island_height)]

     For each island:
         current_shelf = shelves.back()

         If island.width fits in current shelf:
             Place island at (current_shelf.x, current_shelf.y)
             current_shelf.x += island.width + margin
         Else:
             Create new shelf:
             new_y = current_shelf.y + current_shelf.height + margin
             shelves.append((x=0, y=new_y, height=island.height))
             Place island at (0, new_y)

4. Move islands to packed positions:
     For each island:
         offset_x = target_x - current_min_u
         offset_y = target_y - current_min_v

         For each UV in island:
             u += offset_x
             v += offset_y

5. Scale to fit [0,1]²:
     Find max_width, max_height of packed result
     scale = 1.0 / max(max_width, max_height)

     For all UVs:
         u *= scale
         v *= scale
```

### Expected Coverage

- Shelf packing: > 60%
- Optimal packing: > 75% (bonus)

### Bonus: Advanced Packing

For +10 bonus points, implement:
- **MaxRects:** Rectangle packing with multiple heuristics
- **Skyline:** Similar to shelf but with variable heights

---

## Quality Metrics

### Stretch

Measures how much the UV mapping distorts the mesh.

**Computation:**
```
For each triangle:
    Build Jacobian J (2×3 matrix mapping UV → 3D)
    Compute singular values σ1, σ2 via SVD
    stretch = max(σ1/σ2, σ2/σ1)

max_stretch = max over all triangles
```

**Interpretation:**
- Stretch = 1.0: Perfect (no distortion)
- Stretch < 1.3: Excellent
- Stretch < 1.5: Good
- Stretch > 2.0: Poor

### Coverage

Percentage of [0,1]² texture space covered by UVs.

**Computation:**
```
1. Rasterize UVs to 1024×1024 grid
2. Count pixels covered by triangles
3. coverage = pixels_covered / (1024 * 1024)
```

**Interpretation:**
- Coverage > 75%: Excellent
- Coverage > 60%: Good
- Coverage < 50%: Poor (wasting texture space)

---

## References

1. **LSCM Paper:**
   "Least Squares Conformal Maps for Automatic Texture Atlas Generation"
   Lévy et al., SIGGRAPH 2002

2. **Eigen Library:**
   https://eigen.tuxfamily.org/

3. **Mesh Topology:**
   "Polygon Mesh Processing" by Botsch et al.

---

## Tips for Implementation

1. **Start Simple:** Test each component independently
2. **Use Reference Code:** Study `lscm_matrix_example.cpp`
3. **Validate Early:** Use `compare_uvs.py` frequently
4. **Check Memory:** Run valgrind to catch leaks
5. **Profile:** If slow, profile before optimizing

Good luck!
