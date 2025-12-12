"""
Python bindings to C++ UV unwrapping library

This file uses ctypes to wrap the C++ shared library produced in Part 1.
It implements:
 - find_library() to locate the shared library
 - conversion helpers between CMesh and Python Mesh
 - load_mesh(filename)
 - save_mesh(mesh, filename)
 - unwrap(path_or_mesh, params=None)
"""

import ctypes
import os
from pathlib import Path
import numpy as np

# -----------------------------------------------------------------------------
# Find compiled C++ library
# -----------------------------------------------------------------------------
def find_library():
    """
    Find the compiled C++ library and return a Path to it.

    Search order (relative to this file):
      - ../part1_cpp/build/libuvunwrap.so      (Linux default)
      - ../part1_cpp/build/libuvunwrap.dylib    (macOS)
      - ../part1_cpp/build/Release/uvunwrap.dll (Windows/VS)
    """
    import platform
    sys_plat = platform.system().lower()

    names = []
    if sys_plat == "linux":
        names = ["libuvunwrap.so"]
    elif sys_plat == "darwin":
        names = ["libuvunwrap.dylib", "libuvunwrap.so"]
    elif sys_plat.startswith("win"):
        names = ["uvunwrap.dll"]
    else:
        names = ["libuvunwrap.so", "libuvunwrap.dylib", "uvunwrap.dll"]

    here = Path(__file__).resolve().parent

    # candidate relative paths to check
    rel_candidates = [
        Path("..") / "part1_cpp" / "build",
        Path("..") / ".." / "part1_cpp" / "build",
        Path("..") / "part1_cpp" / "build" / "Release",
        Path("..") / ".." / "part1_cpp" / "build" / "Release",
    ]

    for rel in rel_candidates:
        base = (here / rel).resolve()
        for n in names:
            p = base / n
            if p.exists():
                return p

    # also try absolute common locations relative to repo root
    repo_root = here.resolve().parents[2]  # ../../..
    for n in names:
        p = repo_root / "starter_code" / "part1_cpp" / "build" / n
        if p.exists():
            return p

    # If not found, raise a helpful error listing all attempted paths
    tried = []
    for rel in rel_candidates:
        for n in names:
            tried.append(str((here / rel / n).resolve()))
    tried.append(str((repo_root / "starter_code" / "part1_cpp" / names[0]).resolve()))

    raise FileNotFoundError(
        "Could not find libuvunwrap library. Tried:\n  " + "\n  ".join(tried) +
        "\n\nBuild part1_cpp first (cmake .. && make) and ensure the library exists."
    )

# -----------------------------------------------------------------------------
# Define C structures matching C headers
# -----------------------------------------------------------------------------
class CMesh(ctypes.Structure):
    _fields_ = [
        ('vertices', ctypes.POINTER(ctypes.c_float)),
        ('num_vertices', ctypes.c_int),
        ('triangles', ctypes.POINTER(ctypes.c_int)),
        ('num_triangles', ctypes.c_int),
        ('uvs', ctypes.POINTER(ctypes.c_float)),
    ]


class CUnwrapParams(ctypes.Structure):
    _fields_ = [
        ('angle_threshold', ctypes.c_float),
        ('min_island_faces', ctypes.c_int),
        ('pack_islands', ctypes.c_int),
        ('island_margin', ctypes.c_float),
    ]


class CUnwrapResult(ctypes.Structure):
    _fields_ = [
        ('num_islands', ctypes.c_int),
        ('face_island_ids', ctypes.POINTER(ctypes.c_int)),
        ('avg_stretch', ctypes.c_float),
        ('max_stretch', ctypes.c_float),
        ('coverage', ctypes.c_float),
    ]

# -----------------------------------------------------------------------------
# Python Mesh wrapper
# -----------------------------------------------------------------------------
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
        return int(len(self.vertices))

    @property
    def num_triangles(self):
        return int(len(self.triangles))

# -----------------------------------------------------------------------------
# Load library and declare function signatures
# -----------------------------------------------------------------------------
_lib_path = find_library()
_lib = ctypes.CDLL(str(_lib_path))

# Mesh functions
_lib.load_obj.argtypes = [ctypes.c_char_p]
_lib.load_obj.restype = ctypes.POINTER(CMesh)

_lib.save_obj.argtypes = [ctypes.POINTER(CMesh), ctypes.c_char_p]
_lib.save_obj.restype = ctypes.c_int

_lib.free_mesh.argtypes = [ctypes.POINTER(CMesh)]
_lib.free_mesh.restype = None

# Unwrap API
_lib.unwrap_mesh.argtypes = [
    ctypes.POINTER(CMesh),
    ctypes.POINTER(CUnwrapParams),
    ctypes.POINTER(ctypes.POINTER(CUnwrapResult))
]
_lib.unwrap_mesh.restype = ctypes.POINTER(CMesh)

_lib.free_unwrap_result.argtypes = [ctypes.POINTER(CUnwrapResult)]
_lib.free_unwrap_result.restype = None

_lib.compute_quality_metrics.argtypes = [ctypes.POINTER(CMesh), ctypes.POINTER(CUnwrapResult)]
_lib.compute_quality_metrics.restype = None

# -----------------------------------------------------------------------------
# Helper conversions between ctypes CMesh and Python Mesh
# -----------------------------------------------------------------------------
def _cmesh_to_python(cm_ptr):
    """Convert ctypes POINTER(CMesh) to Python Mesh object (copies data)."""
    if not cm_ptr:
        return None
    cm = cm_ptr.contents
    V = int(cm.num_vertices)
    F = int(cm.num_triangles)

    # vertices (flat 3*V)
    if V > 0:
        verts_arr = np.ctypeslib.as_array(cm.vertices, shape=(3 * V,))
        verts = verts_arr.reshape((V, 3)).copy()
    else:
        verts = np.zeros((0, 3), dtype=np.float32)

    # triangles (flat 3*F)
    if F > 0:
        tris_arr = np.ctypeslib.as_array(cm.triangles, shape=(3 * F,))
        tris = tris_arr.reshape((F, 3)).copy()
    else:
        tris = np.zeros((0, 3), dtype=np.int32)

    # uvs may be NULL
    uvs = None
    if bool(cm.uvs):
        uvs_arr = np.ctypeslib.as_array(cm.uvs, shape=(2 * V,))
        uvs = uvs_arr.reshape((V, 2)).copy()

    return Mesh(verts, tris, uvs)


def _python_to_cmesh(py_mesh):
    """Create a CMesh instance (ctypes) from Python Mesh.
    Returns (cm_struct, buffers) where buffers list must be kept alive while used.
    """
    V = py_mesh.num_vertices
    F = py_mesh.num_triangles

    # create ctypes arrays (these own the memory)
    vert_flat = py_mesh.vertices.reshape(-1).astype(np.float32)
    tri_flat = py_mesh.triangles.reshape(-1).astype(np.int32)

    vert_buf = (ctypes.c_float * (3 * V))(*vert_flat)
    tri_buf = (ctypes.c_int * (3 * F))(*tri_flat)

    if py_mesh.uvs is not None:
        uv_flat = py_mesh.uvs.reshape(-1).astype(np.float32)
        uv_buf = (ctypes.c_float * (2 * V))(*uv_flat)
        uv_ptr = ctypes.cast(uv_buf, ctypes.POINTER(ctypes.c_float))
    else:
        uv_buf = None
        uv_ptr = ctypes.POINTER(ctypes.c_float)()

    cm = CMesh()
    cm.vertices = ctypes.cast(vert_buf, ctypes.POINTER(ctypes.c_float))
    cm.num_vertices = ctypes.c_int(V)
    cm.triangles = ctypes.cast(tri_buf, ctypes.POINTER(ctypes.c_int))
    cm.num_triangles = ctypes.c_int(F)
    cm.uvs = uv_ptr

    # Return struct and buffers to keep alive (caller must retain buffers)
    return cm, (vert_buf, tri_buf, uv_buf)

# -----------------------------------------------------------------------------
# Public wrappers
# -----------------------------------------------------------------------------
def load_mesh(filename):
    """
    Load mesh from OBJ file using the C library and return Python Mesh (copies data).
    """
    p = Path(filename)
    if not p.exists():
        raise FileNotFoundError(filename)

    cptr = _lib.load_obj(str(p).encode('utf-8'))
    if not cptr:
        raise RuntimeError(f"load_obj failed for {filename}")

    py = _cmesh_to_python(cptr)
    _lib.free_mesh(cptr)
    return py


def save_mesh(mesh, filename):
    """
    Save Python Mesh to OBJ via C library.
    """
    # convert to CMesh and call save_obj
    cm, buffers = _python_to_cmesh(mesh)
    # call save_obj with pointer to cm
    res = _lib.save_obj(ctypes.byref(cm), str(filename).encode('utf-8'))
    if res != 0:
        raise RuntimeError(f"save_obj failed (code {res})")
    return res


def unwrap(path_or_mesh, params=None):
    """
    Unwrap a mesh using the C++ library.

    Accepts either a file path (string/Path) or a Python Mesh object.
    Returns (py_mesh, metrics_dict).
    """
    # default params
    default = {
        'angle_threshold': 30.0,
        'min_island_faces': 5,
        'pack_islands': True,
        'island_margin': 0.02
    }
    if params:
        default.update(params)

    # prepare UnwrapParams C struct
    cparams = CUnwrapParams()
    cparams.angle_threshold = ctypes.c_float(default['angle_threshold'])
    cparams.min_island_faces = ctypes.c_int(default['min_island_faces'])
    cparams.pack_islands = ctypes.c_int(1 if default['pack_islands'] else 0)
    cparams.island_margin = ctypes.c_float(default['island_margin'])

    need_free_input = False

    # Load mesh into C (if path) or convert Python mesh to C
    if isinstance(path_or_mesh, (str, Path)):
        p = Path(path_or_mesh)
        if not p.exists():
            raise FileNotFoundError(path_or_mesh)
        c_in = _lib.load_obj(str(p).encode('utf-8'))
        if not c_in:
            raise RuntimeError(f"load_obj failed for {path_or_mesh}")
        need_free_input = True
    elif isinstance(path_or_mesh, Mesh):
        # convert Python mesh to CMesh and keep buffers alive
        cm, buffers = _python_to_cmesh(path_or_mesh)
        # create pointer to cm
        c_in = ctypes.pointer(cm)
        # Keep buffers alive by attaching to pointer (so GC doesn't free)
        c_in._buffers = buffers
    else:
        raise TypeError("path_or_mesh must be a file path or Mesh instance")

    # prepare pointer to receive result pointer
    result_ptr = ctypes.POINTER(CUnwrapResult)()
    out_mesh_ptr = _lib.unwrap_mesh(c_in, ctypes.byref(cparams), ctypes.byref(result_ptr))
    if not out_mesh_ptr:
        if need_free_input:
            _lib.free_mesh(c_in)
        raise RuntimeError("unwrap_mesh failed")

    # Convert returned mesh to Python
    py_out = _cmesh_to_python(out_mesh_ptr)

    # Convert result struct to dict
    res = {}
    if bool(result_ptr):
        r = result_ptr.contents
        res = {
            'num_islands': int(r.num_islands),
            'avg_stretch': float(r.avg_stretch),
            'max_stretch': float(r.max_stretch),
            'coverage': float(r.coverage)
        }

    # Free C resources
    try:
        if bool(result_ptr):
            _lib.free_unwrap_result(result_ptr)
    except Exception:
        # ignore freeing errors but log nothing here (caller handles)
        pass

    _lib.free_mesh(out_mesh_ptr)
    if need_free_input:
        _lib.free_mesh(c_in)

    return py_out, res

# -----------------------------------------------------------------------------
# Quick test when run as script
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    print("Bindings test: library at", find_library())
    # quick smoke test if test mesh exists
    test_mesh = Path(__file__).resolve().parents[2] / "test_data" / "meshes" / "01_cube.obj"
    if test_mesh.exists():
        try:
            m, metrics = unwrap(str(test_mesh))
            print("Unwrapped:", m.num_vertices, "verts", m.num_triangles, "tris")
            print("Metrics:", metrics)
        except Exception as e:
            print("Unwrap failed:", e)
    else:
        print("No test mesh found at", test_mesh)

