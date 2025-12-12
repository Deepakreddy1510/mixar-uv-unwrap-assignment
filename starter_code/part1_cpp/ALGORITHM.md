# ALGORITHM.md

## Overview
This implementation provides an automatic UV unwrapping pipeline consisting of:
- Topology builder: extract unique edges and adjacent faces (edge_face adjacency).
- Seam detection: hybrid strategy using dihedral angles, boundary-loop representatives,
  and a fallback. Also uses angular defect heuristics for curvature-aware refinement.
- Island extraction: connected components in face graph after seam cuts.
- Parameterization: LSCM per-island (using Eigen for linear algebra).
- Packing: simple shelf-like packing with margins and scaling to [0,1]².
- Quality metrics: stretch (SVD-based), coverage (rasterized), and angle distortion.

## Topology
- Build edges for every triangle (normalize edge as (min,max)).
- Maintain `edge_faces` array: for each edge store up to two adjacent faces; boundary edge has -1.
- Validate: compute Euler characteristic `χ = V - E + F` for sanity checks.

## Seam detection
- For each internal edge compute dihedral angle between adjacent face normals (degrees).
- Mark edges whose dihedral ≥ threshold (adaptive thresholding included).
- Add boundary-edge representatives: group boundary edges by connected components (loop) and pick one per loop.
- Fallback: if no seams, pick edge with maximal dihedral.
- Output: list of seam edge indices.

## Island extraction
- Remove seam edges and build face adjacency graph using remaining edges.
- Compute connected components (BFS) over faces.
- Assign island id per face, return `face_island_ids` and `num_islands`.

## LSCM Parameterization
- For each island:
  - Build local vertex index mapping (global → local).
  - Build local triangle list and call `lscm_parameterize`.
  - Copy per-vertex UVs into global result (respecting global indices).
- Skip islands smaller than `min_island_faces`.
- Use Eigen sparse linear solver for LSCM.

## Packing
- Compute per-island UV bounding boxes.
- Sort islands by height and place them into rows (shelf algorithm).
- Apply uniform scale to fit all islands into [0,1]² while respecting margin.

## Quality Metrics
- **Stretch**: per-triangle Jacobian (3×2) singular values via SVD, report max ratio σ1/σ2.
- **Coverage**: rasterize triangles into resolution×resolution grid and report filled pixel fraction.
- **Angle distortion**: per-triangle compare 3D angles vs UV angles, report max absolute difference.

## Notes / Limitations
- Islands that are very small are skipped for LSCM (left as zero-uvs or rely on packing).
- Packing is shelf-based and not optimal; fine for baseline functionality.
- Provided heuristics for seam detection work for test meshes; may need tuning for complex geometry.


