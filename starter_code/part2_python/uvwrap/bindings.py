"""
Python bindings to C++ UV unwrapping library

TEMPLATE - YOU IMPLEMENT

Uses ctypes to wrap the C++ shared library.
Alternative: Use pybind11 for cleaner bindings (bonus points).
"""

import ctypes
import os
from pathlib import Path
import numpy as np

# TODO: Find the compiled shared library
# Look in: ../part1_cpp/build/libuvunwrap.so (Linux)
#          ../part1_cpp/build/libuvunwrap.dylib (macOS)
#          ../part1_cpp/build/Release/uvunwrap.dll (Windows)

def find_library():
    """
    Find the compiled C++ library

    Returns:
        Path to library file
    """
    # TODO: Implement library finding logic
    #
    # Hints:
    # 1. Look relative to this file's location
    # 2. Check common build directories
    # 3. Handle different platforms (Linux, macOS, Windows)
    # 4. Raise error if not found

    pass  # YOUR CODE HERE


# Load library
# TODO: Uncomment when implemented
# _lib_path = find_library()
# _lib = ctypes.CDLL(str(_lib_path))


# Define C structures matching mesh.h
class CMesh(ctypes.Structure):
    """
    Matches the Mesh struct in mesh.h
    """
    _fields_ = [
        ('vertices', ctypes.POINTER(ctypes.c_float)),
        ('num_vertices', ctypes.c_int),
        ('triangles', ctypes.POINTER(ctypes.c_int)),
        ('num_triangles', ctypes.c_int),
        ('uvs', ctypes.POINTER(ctypes.c_float)),
    ]


class CUnwrapParams(ctypes.Structure):
    """
    Matches UnwrapParams struct in unwrap.h
    """
    _fields_ = [
        ('angle_threshold', ctypes.c_float),
        ('min_island_faces', ctypes.c_int),
        ('pack_islands', ctypes.c_int),
        ('island_margin', ctypes.c_float),
    ]


class CUnwrapResult(ctypes.Structure):
    """
    Matches UnwrapResult struct in unwrap.h
    """
    _fields_ = [
        ('num_islands', ctypes.c_int),
        ('face_island_ids', ctypes.POINTER(ctypes.c_int)),
        ('avg_stretch', ctypes.c_float),
        ('max_stretch', ctypes.c_float),
        ('coverage', ctypes.c_float),
    ]


# TODO: Define function signatures
#
# Example:
# _lib.load_obj.argtypes = [ctypes.c_char_p]
# _lib.load_obj.restype = ctypes.POINTER(CMesh)
#
# _lib.free_mesh.argtypes = [ctypes.POINTER(CMesh)]
# _lib.free_mesh.restype = None
#
# ... etc for all functions


class Mesh:
    """
    Python wrapper for C mesh

    Attributes:
        vertices: numpy array (N, 3) of vertex positions
        triangles: numpy array (M, 3) of triangle indices
        uvs: numpy array (N, 2) of UV coordinates (optional)
    """

    def __init__(self, vertices, triangles, uvs=None):
        self.vertices = np.array(vertices, dtype=np.float32)
        self.triangles = np.array(triangles, dtype=np.int32)
        self.uvs = np.array(uvs, dtype=np.float32) if uvs is not None else None

    @property
    def num_vertices(self):
        return len(self.vertices)

    @property
    def num_triangles(self):
        return len(self.triangles)


def load_mesh(filename):
    """
    Load mesh from OBJ file

    Args:
        filename: Path to OBJ file

    Returns:
        Mesh object

    IMPLEMENTATION REQUIRED
    """
    # TODO: Implement
    #
    # Steps:
    # 1. Call C library load_obj function
    # 2. Convert C mesh to Python Mesh object
    # 3. Copy data from C arrays to numpy arrays
    # 4. Free C mesh
    # 5. Return Python Mesh

    pass  # YOUR CODE HERE


def save_mesh(mesh, filename):
    """
    Save mesh to OBJ file

    Args:
        mesh: Mesh object
        filename: Output path

    IMPLEMENTATION REQUIRED
    """
    # TODO: Implement
    #
    # Steps:
    # 1. Convert Python Mesh to C mesh
    # 2. Call C library save_obj function
    # 3. Free C mesh

    pass  # YOUR CODE HERE


def unwrap(mesh, params=None):
    """
    Unwrap mesh using LSCM

    Args:
        mesh: Mesh object
        params: Dictionary of parameters:
            - angle_threshold: float (default 30.0)
            - min_island_faces: int (default 10)
            - pack_islands: bool (default True)
            - island_margin: float (default 0.02)

    Returns:
        tuple: (unwrapped_mesh, result_dict)
            unwrapped_mesh: Mesh with UVs
            result_dict: {
                'num_islands': int,
                'max_stretch': float,
                'avg_stretch': float,
                'coverage': float,
            }

    IMPLEMENTATION REQUIRED
    """
    # TODO: Implement
    #
    # Steps:
    # 1. Convert Python Mesh to C mesh
    # 2. Set up C parameters
    # 3. Call C library unwrap_mesh function
    # 4. Convert result C mesh to Python Mesh
    # 5. Extract result metrics
    # 6. Free C meshes
    # 7. Return Python objects

    pass  # YOUR CODE HERE


# Example usage (for testing)
if __name__ == "__main__":
    # Test loading
    print("Testing bindings...")

    # TODO: Test with a simple mesh
    # mesh = load_mesh("../../test_data/meshes/01_cube.obj")
    # print(f"Loaded: {mesh.num_vertices} vertices, {mesh.num_triangles} triangles")

    # result_mesh, metrics = unwrap(mesh)
    # print(f"Unwrapped: {metrics['num_islands']} islands")

    pass
