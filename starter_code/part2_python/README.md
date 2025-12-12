# Part 2 — Python Batch Processor

## Overview
This module provides Python bindings to the C++ unwrapping library and a CLI tool:
- `uvwrap/bindings.py` — ctypes bindings: `load_mesh`, `save_mesh`, `unwrap`.
- `uvwrap/metrics.py` — quality metrics: `compute_stretch`, `compute_coverage`, `compute_angle_distortion`.
- `cli.py` — command-line interface:
  - `unwrap` — unwrap a single mesh and save output.
  - `batch` — multi-threaded batch processing of a folder of OBJ files.
  - `optimize` — simple grid-search parameter optimization.
  - `analyze` — compute metrics and print JSON report.

## Quick start
1. Build Part 1 (C++):
   ```bash
   cd starter_code/part1_cpp
   mkdir -p build && cd build
   cmake ..
   make -j$(nproc)

