"""
Quality metrics for UV mappings

TEMPLATE - YOU IMPLEMENT

Computes:
- Stretch: UV→3D Jacobian singular value ratio
- Coverage: Percentage of [0,1]² used
- Angle distortion: Max angle difference

See reference/metrics_spec.md for exact formulas
"""

import numpy as np


def compute_stretch(mesh, uvs):
    """
    Compute maximum stretch across all triangles

    Stretch measures how much the UV mapping distorts the mesh.

    Algorithm:
    For each triangle:
        1. Build Jacobian J (UV → 3D mapping)
        2. Compute singular values σ1, σ2 via SVD
        3. stretch = max(σ1/σ2, σ2/σ1)

    Args:
        mesh: Mesh object with vertices and triangles
        uvs: UV coordinates (N, 2) array

    Returns:
        float: Maximum stretch value

    IMPLEMENTATION REQUIRED
    See reference/metrics_spec.md for exact formula
    """
    # TODO: Implement
    #
    # Steps:
    # 1. For each triangle:
    #      Get 3D positions p0, p1, p2
    #      Get UV positions uv0, uv1, uv2
    #
    # 2. Compute edge vectors:
    #      dp1 = p1 - p0 (3D)
    #      dp2 = p2 - p0 (3D)
    #      duv1 = uv1 - uv0 (2D)
    #      duv2 = uv2 - uv0 (2D)
    #
    # 3. Build Jacobian J:
    #      J = [dp1 | dp2] @ inv([[duv1], [duv2]])
    #      Result is 3×2 matrix
    #
    # 4. Compute SVD:
    #      U, S, Vt = np.linalg.svd(J)
    #      σ1, σ2 = S[0], S[1]
    #
    # 5. Stretch for this triangle:
    #      stretch = max(σ1/σ2, σ2/σ1)
    #
    # 6. Return max stretch across all triangles

    pass  # YOUR CODE HERE


def compute_coverage(uvs, triangles, resolution=1024):
    """
    Compute UV coverage (percentage of [0,1]² covered)

    Algorithm:
    1. Rasterize UVs to resolution×resolution grid
    2. Mark pixels covered by triangles
    3. Return percentage of pixels covered

    Args:
        uvs: UV coordinates (N, 2) array
        triangles: Triangle indices (M, 3) array
        resolution: Grid resolution (default 1024)

    Returns:
        float: Coverage (0.0 to 1.0)

    IMPLEMENTATION REQUIRED
    See reference/metrics_spec.md for details
    """
    # TODO: Implement
    #
    # Steps:
    # 1. Create boolean grid (resolution × resolution)
    # 2. For each triangle:
    #      Scale UVs to grid coordinates
    #      Rasterize triangle (mark covered pixels)
    # 3. Count covered pixels
    # 4. Return covered / total

    # Helper function for rasterization
    def rasterize_triangle(grid, uv0, uv1, uv2):
        """Mark pixels covered by triangle"""
        # TODO: Implement triangle rasterization
        #
        # Hints:
        # - Find bounding box
        # - For each pixel in bounding box:
        #     Check if inside triangle (barycentric test)
        #     If inside, mark grid[y, x] = True
        pass

    pass  # YOUR CODE HERE


def compute_angle_distortion(mesh, uvs):
    """
    Compute maximum angle distortion

    Measures how much angles change between 3D and UV space.

    Algorithm:
    For each triangle:
        Compute 3 angles in 3D
        Compute 3 angles in UV
        Find max |angle_3d - angle_uv|

    Args:
        mesh: Mesh object
        uvs: UV coordinates (N, 2) array

    Returns:
        float: Max angle distortion (radians)

    IMPLEMENTATION REQUIRED
    See reference/metrics_spec.md for formula
    """
    # TODO: Implement
    #
    # Steps:
    # 1. For each triangle:
    #      Compute 3 angles in 3D (using dot product)
    #      Compute 3 angles in UV (using dot product)
    #      Find max difference
    # 2. Return max across all triangles

    def compute_triangle_angles_3d(p0, p1, p2):
        """Compute 3 angles of 3D triangle"""
        # TODO: Implement
        # Angle at p0: arccos(dot(normalize(p1-p0), normalize(p2-p0)))
        # Similar for p1 and p2
        pass

    def compute_triangle_angles_2d(uv0, uv1, uv2):
        """Compute 3 angles of 2D triangle"""
        # TODO: Implement (same as 3D but in 2D)
        pass

    pass  # YOUR CODE HERE


# Example usage
if __name__ == "__main__":
    # Test with simple triangle
    vertices = np.array([
        [0, 0, 0],
        [1, 0, 0],
        [0.5, 0.866, 0],
    ], dtype=np.float32)

    triangles = np.array([[0, 1, 2]], dtype=np.int32)

    uvs = np.array([
        [0, 0],
        [1, 0],
        [0.5, 0.866],
    ], dtype=np.float32)

    class SimpleMesh:
        def __init__(self, vertices, triangles):
            self.vertices = vertices
            self.triangles = triangles

    mesh = SimpleMesh(vertices, triangles)

    # Should be 1.0 (no distortion for this triangle)
    print(f"Stretch: {compute_stretch(mesh, uvs)}")
    print(f"Coverage: {compute_coverage(uvs, triangles)}")
    print(f"Angle distortion: {compute_angle_distortion(mesh, uvs)}")
