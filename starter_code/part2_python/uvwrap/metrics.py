"""
Quality metrics for UV mappings

Implements:
- compute_stretch(mesh, uvs)
- compute_coverage(uvs, triangles, resolution=1024)
- compute_angle_distortion(mesh, uvs)

Notes:
- Mesh is expected to have attributes:
    .vertices -> (N,3) numpy array
    .triangles -> (M,3) numpy array (or flat length 3*M)
- UVs is (N,2) numpy array
"""

import numpy as np


def _ensure_tris_array(triangles):
    arr = np.array(triangles, copy=False)
    if arr.ndim == 1 and arr.size % 3 == 0:
        arr = arr.reshape((-1, 3))
    return arr.astype(np.int32)


def compute_stretch(mesh, uvs):
    """
    Compute maximum stretch across all triangles.

    Return: float (max stretch across triangles).
    """
    verts = np.array(mesh.vertices, dtype=np.float64)
    tris = _ensure_tris_array(mesh.triangles)

    uvs = np.array(uvs, dtype=np.float64)

    max_stretch = 0.0

    for tri in tris:
        i0, i1, i2 = tri
        p0 = verts[i0]; p1 = verts[i1]; p2 = verts[i2]
        uv0 = uvs[i0]; uv1 = uvs[i1]; uv2 = uvs[i2]

        # Edge vectors in 3D
        dp1 = p1 - p0  # 3
        dp2 = p2 - p0  # 3

        # Edge vectors in UV (2D)
        duv1 = uv1 - uv0  # 2
        duv2 = uv2 - uv0  # 2

        # Build 3x2 and 2x2 matrices
        A = np.column_stack((dp1, dp2))  # 3x2
        B = np.column_stack((duv1, duv2))  # 2x2

        # If UV triangle is degenerate, skip (cannot invert)
        det = B[0, 0] * B[1, 1] - B[0, 1] * B[1, 0]
        if abs(det) < 1e-12:
            # degenerate UV mapping — treat as high stretch
            continue

        try:
            Binv = np.linalg.inv(B)  # 2x2
        except np.linalg.LinAlgError:
            continue

        # Jacobian J: 3x2 = A (3x2) @ Binv (2x2)
        J = A.dot(Binv)

        # Compute singular values of J (only two non-zero singular values at most)
        # Use SVD; S contains singular values sorted descending
        try:
            s = np.linalg.svd(J, compute_uv=False)
        except np.linalg.LinAlgError:
            continue

        if s.size < 2:
            continue
        s1, s2 = s[0], s[1]
        if s2 <= 0:
            continue
        stretch = max(s1 / s2, s2 / s1)
        if stretch > max_stretch:
            max_stretch = stretch

    # if no triangles produced a valid value, return 1.0 (no distortion)
    if max_stretch == 0.0:
        return 1.0
    return float(max_stretch)


def compute_coverage(uvs, triangles, resolution=1024):
    """
    Rasterize UV triangles into resolution x resolution grid and compute coverage fraction.

    Returns: float in [0,1]
    """
    uvs = np.array(uvs, dtype=np.float64)
    tris = _ensure_tris_array(triangles)

    if uvs.size == 0 or tris.size == 0:
        return 0.0

    # Normalize UVs to [0,1]
    # Some UVs may already be in range; clamp afterwards.
    uv = np.array(uvs, dtype=np.float64)
    # If UVs are outside 0..1, clamp them to that range for rasterization
    uv = np.clip(uv, 0.0, 1.0)

    # grid boolean
    H = W = int(resolution)
    grid = np.zeros((H, W), dtype=np.bool_)

    def rasterize_triangle(grid, uv0, uv1, uv2):
        # Scale to grid coordinates [0, W-1] and [0, H-1], using pixel centers
        pts = np.array([uv0, uv1, uv2], dtype=np.float64)
        pts_pix = pts * np.array([W - 1.0, H - 1.0])
        xs = pts_pix[:, 0]
        ys = pts_pix[:, 1]

        # bounding box
        min_x = int(max(0, np.floor(xs.min())))
        max_x = int(min(W - 1, np.ceil(xs.max())))
        min_y = int(max(0, np.floor(ys.min())))
        max_y = int(min(H - 1, np.ceil(ys.max())))
        if min_x > max_x or min_y > max_y:
            return

        # Precompute for barycentric
        x0, y0 = xs[0], ys[0]
        x1, y1 = xs[1], ys[1]
        x2, y2 = xs[2], ys[2]
        denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        if abs(denom) < 1e-12:
            return  # degenerate triangle in pixel space

        # iterate pixels in bbox
        for py in range(min_y, max_y + 1):
            for px in range(min_x, max_x + 1):
                # pixel center
                cx = px + 0.5
                cy = py + 0.5
                # barycentric coordinates
                w0 = ((y1 - y2) * (cx - x2) + (x2 - x1) * (cy - y2)) / denom
                w1 = ((y2 - y0) * (cx - x2) + (x0 - x2) * (cy - y2)) / denom
                w2 = 1.0 - w0 - w1
                # inside if all barycentric in [0,1]
                if (w0 >= -1e-8) and (w1 >= -1e-8) and (w2 >= -1e-8):
                    grid[py, px] = True

    # Rasterize each triangle
    for tri in tris:
        i0, i1, i2 = tri
        uv0 = uv[i0]; uv1 = uv[i1]; uv2 = uv[i2]
        # If triangle is completely outside [0,1], skip
        if ((uv0[0] < 0 and uv1[0] < 0 and uv2[0] < 0) or
            (uv0[0] > 1 and uv1[0] > 1 and uv2[0] > 1) or
            (uv0[1] < 0 and uv1[1] < 0 and uv2[1] < 0) or
            (uv0[1] > 1 and uv1[1] > 1 and uv2[1] > 1)):
            continue
        rasterize_triangle(grid, uv0, uv1, uv2)

    covered = np.count_nonzero(grid)
    total = float(H * W)
    return float(covered) / total if total > 0 else 0.0


def compute_angle_distortion(mesh, uvs):
    """
    Compute maximum angle distortion (in radians) between 3D triangle angles and UV triangle angles.
    """
    verts = np.array(mesh.vertices, dtype=np.float64)
    tris = _ensure_tris_array(mesh.triangles)
    uvs = np.array(uvs, dtype=np.float64)

    def tri_angles_3d(p0, p1, p2):
        # return three angles at p0,p1,p2
        a0 = p1 - p0
        b0 = p2 - p0
        a1 = p0 - p1
        b1 = p2 - p1
        a2 = p0 - p2
        b2 = p1 - p2

        def ang(u, v):
            nu = np.linalg.norm(u)
            nv = np.linalg.norm(v)
            if nu < 1e-12 or nv < 1e-12:
                return 0.0
            cosv = np.dot(u, v) / (nu * nv)
            cosv = max(-1.0, min(1.0, cosv))
            return np.arccos(cosv)

        return ang(a0, b0), ang(a1, b1), ang(a2, b2)

    def tri_angles_2d(u0, u1, u2):
        a0 = u1 - u0
        b0 = u2 - u0
        a1 = u0 - u1
        b1 = u2 - u1
        a2 = u0 - u2
        b2 = u1 - u2

        def ang(u, v):
            nu = np.hypot(u[0], u[1])
            nv = np.hypot(v[0], v[1])
            if nu < 1e-12 or nv < 1e-12:
                return 0.0
            cosv = (u[0]*v[0] + u[1]*v[1]) / (nu * nv)
            cosv = max(-1.0, min(1.0, cosv))
            return np.arccos(cosv)

        return ang(a0, b0), ang(a1, b1), ang(a2, b2)

    max_diff = 0.0
    for tri in tris:
        i0, i1, i2 = tri
        p0 = verts[i0]; p1 = verts[i1]; p2 = verts[i2]
        u0 = uvs[i0]; u1 = uvs[i1]; u2 = uvs[i2]

        angles3 = tri_angles_3d(p0, p1, p2)
        angles2 = tri_angles_2d(u0, u1, u2)

        # compute max absolute difference across the three angles
        diffs = [abs(a3 - a2) for a3, a2 in zip(angles3, angles2)]
        local_max = max(diffs)
        if local_max > max_diff:
            max_diff = local_max

    return float(max_diff)


# Example usage (quick smoke test)
if __name__ == "__main__":
    vertices = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0.5, 0.8660254, 0],
    ], dtype=np.float32)
    triangles = np.array([[0, 1, 2]], dtype=np.int32)
    uvs = np.array([
        [0, 0],
        [1, 0],
        [0.5, 0.8660254],
    ], dtype=np.float32)

    class SimpleMesh:
        def __init__(self, vertices, triangles):
            self.vertices = vertices
            self.triangles = triangles

    mesh = SimpleMesh(vertices, triangles)
    print("Stretch:", compute_stretch(mesh, uvs))
    print("Coverage:", compute_coverage(uvs, triangles, resolution=128))
    print("Angle distortion (rad):", compute_angle_distortion(mesh, uvs))

