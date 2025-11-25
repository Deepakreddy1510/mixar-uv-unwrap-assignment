# Mixar SDE Technical Assignment

**Project:** Automatic UV Unwrapping with Seam Detection

**Time Estimate:** 18-23 hours

**Deadline:** 7 days from receipt

---

## Overview

Build a production-ready automatic UV unwrapping system that segments 3D meshes into UV islands, detects optimal seams, and performs LSCM (Least Squares Conformal Maps) parameterization.

### What You'll Build

- **Part 1 (C++):** Core unwrapping engine with topology analysis, seam detection, and LSCM parameterization (70 points)
- **Part 2 (Python):** Multi-threaded batch processor with quality metrics and parameter optimization (35 points)
- **Part 3 (Blender):** Production addon with caching, seam editing, and live preview (35 points)

**Total: 140 points** (Pass threshold: 100/140 = 71%)

---

## Part 1: C++ Unwrapping Engine (70 points)

### Time Estimate
8-10 hours

### What You'll Implement

**Step 1: Topology Builder (1-2 hours)**
- Build edge list from triangle mesh
- Find adjacent faces for each edge
- Validate using Euler characteristic

**Step 2: Seam Detection (2-3 hours)**
- Spanning tree on dual graph
- Angular defect computation
- Seam refinement based on curvature

**Step 3: LSCM Parameterization (4-5 hours)** ← Most complex
- Build LSCM sparse matrix
- Set boundary conditions
- Solve sparse linear system
- Normalize UVs to [0,1]²

**Step 4: Island Packing (1 hour)**
- Compute bounding boxes
- Shelf packing algorithm
- Scale to fit [0,1]²

### Resources Provided

- Complete API specifications in header files
- Reference implementation of LSCM matrix assembly
- Detailed algorithm descriptions in `reference/algorithms.md`
- Mathematical background in `reference/lscm_math.md`
- Complete test suite
- Validation tools

### Expected Results

| Test | Criteria |
|------|----------|
| Topology | Correct edge count, adjacency, Euler validation |
| Seam Detection | Cube: 7-9 seams, Sphere: 1-3 seams |
| LSCM | Max stretch < 1.5, Cylinder < 1.2 |
| Packing | No overlaps, fits in [0,1]², >60% coverage |
| Performance | 10k vertices < 5s |

### Deliverables

```
part1_cpp/
├── src/*.cpp                 # YOUR IMPLEMENTATIONS
├── ALGORITHM.md              # Explain your approach
└── TEST_RESULTS.txt          # Test output showing 10/10 pass
```

---

## Part 2: Python Batch Processor (35 points)

### Time Estimate
4-5 hours

### What You'll Implement

1. **Python Bindings (1 hour)**
   - Wrap C++ library using ctypes or pybind11
   - Handle memory management

2. **Quality Metrics (1.5 hours)**
   - Stretch: UV→3D Jacobian singular value ratio
   - Coverage: Rasterize UVs, compute % of [0,1]² used
   - Angle Distortion: Max angle difference between 3D and UV

3. **Multi-threaded Processor (1 hour)**
   - ThreadPoolExecutor for parallel processing
   - Thread-safe progress tracking

4. **CLI Tool (1 hour)**
   - Commands: unwrap, batch, optimize, analyze
   - Progress bars, formatted output

5. **Parameter Optimization (0.5 hours)**
   - Grid search over parameter space

### Resources Provided

- Exact metric formulas in `reference/metrics_spec.md`
- Complete example implementations
- Test suite

### CLI Commands

```bash
# Single file
python cli.py unwrap input.obj output.obj --angle-threshold 30

# Batch processing
python cli.py batch input_dir/ output_dir/ --threads 8 --report metrics.json

# Optimize parameters
python cli.py optimize input.obj --metric stretch --save-params best.json

# Analyze quality
python cli.py analyze mesh.obj
```

### Deliverables

```
part2_python/
├── uvwrap/*.py               # YOUR IMPLEMENTATIONS
├── cli.py                    # Command-line tool
└── README.md                 # Usage examples
```

---

## Part 3: Blender Add-on (35 points)

### Time Estimate
6-8 hours

### What You'll Implement

1. **Main Unwrap Operator (2 hours)**
   - Extract mesh data directly (no import/export)
   - Call C++ library via Python bindings
   - Apply UVs back to Blender mesh

2. **Caching System (1.5 hours)**
   - Hash mesh geometry + parameters
   - Store results with timestamp
   - Clean up old cache entries

3. **UI Panel (1 hour)**
   - Settings controls
   - Results display
   - Cache status indicator

4. **Seam Editing (1 hour)**
   - Mark/clear seams
   - Invalidate cache when changed

5. **Batch Processing (0.5 hours)**
   - Unwrap all meshes in scene
   - Progress bar

6. **Live Preview (2 hours)**
   - Modal operator with timer
   - Debouncing (200ms)
   - FPS display

### Resources Provided

- UI mockup and design specifications
- Complete implementation guide
- Direct mesh extraction examples
- Caching implementation example

### Deliverables

```
part3_blender/
├── *.py                      # YOUR IMPLEMENTATIONS
├── demo.blend                # Example file
├── screenshots/              # UI screenshots
└── README.md                 # Installation & usage
```

---

## Grading Summary

| Part | Points | Key Components |
|------|--------|----------------|
| Part 1: C++ | 70 | Topology (10), Seams (15), LSCM (20), Packing (5), Performance (5), Quality (10), Docs (5) |
| Part 2: Python | 35 | Bindings (5), Multi-thread (8), Metrics (8), CLI (5), Optimization (6), Errors (3) |
| Part 3: Blender | 35 | Unwrap (5), UI (5), Cache (8), Seams (5), Batch (3), Undo (3), Live (6) |
| **TOTAL** | **140** | |

**Grade Boundaries:**
- 120-140: Excellent
- 100-119: Good
- 90-99: Acceptable
- <90: Does not meet requirements

---

## Submission

**Submit via the provided Google Form link.**

**You will need to provide:**
1. GitHub repository link (private repo, invite: **rahul237** and **Raj-Mixar**)
2. Difficulties encountered / Highlights of your solution in README.md

---

## Getting Started

1. **Fork the repository:**
   ```bash
   # Fork https://github.com/Mixar-AI/mixar-uv-unwrap-assignment
   git clone https://github.com/YOUR_USERNAME/mixar-uv-unwrap-assignment.git
   cd mixar-uv-unwrap-assignment
   ```

2. **Read README.txt** for detailed setup instructions

3. **Build Part 1:**
   ```bash
   cd starter_code/part1_cpp
   mkdir build && cd build
   cmake ..
   make
   ./test_unwrap
   ```

4. **Setup Part 2:**
   ```bash
   cd starter_code/part2_python
   pip install -r requirements.txt
   ```

5. **Test Part 3:**
   - Open Blender (4.2+)
   - Install add-on from `starter_code/part3_blender/`

---

## Resources in Package

### Documentation
- `starter_code/part1_cpp/reference/algorithms.md` - Algorithm descriptions
- `starter_code/part1_cpp/reference/lscm_math.md` - Mathematical derivations
- `starter_code/part2_python/reference/metrics_spec.md` - Exact metric formulas

### Reference Code
- `starter_code/part1_cpp/reference/lscm_matrix_example.cpp` - LSCM matrix assembly

### Test Data
- `test_data/meshes/` - 4 test meshes (cube, cylinder, sphere, torus)

---

## Tips for Success

✅ Start with Part 1 - It's the foundation
✅ Use reference implementations - We provide examples for a reason
✅ Test incrementally - Don't wait until the end
✅ Run the test suite frequently - Catch issues early
✅ Document as you go - Write ALGORITHM.md while coding
✅ Make small commits — Track incremental progress 
✅ Ask questions early - Better than wrong assumptions

---

## Common Pitfalls to Avoid

❌ Don't implement your own sparse matrix solver (use Eigen)
❌ Don't use import/export for Blender (extract data directly)
❌ Don't skip testing (run test_unwrap frequently)
❌ Don't forget memory management (check with valgrind)
❌ Don't ignore edge cases (degenerate triangles, closed islands)

---

## Questions?

**Email:** ajay@mixar.app
**Response time:** Within 1 business day (weekdays)
**Ask within:** First 72 hours for best response

---

Good luck! This assignment reflects real work at Mixar. We're excited to see your solution.

— The Mixar Team
