================================================================================
MIXAR SDE TECHNICAL ASSIGNMENT
UV Unwrapping System
================================================================================

GETTING STARTED
---------------

1. FORK THE ASSIGNMENT REPOSITORY
   https://github.com/Mixar-AI/mixar-uv-unwrap-assignment

2. READ THE ASSIGNMENT
   - Open ASSIGNMENT.md for complete requirements
   - This README provides setup instructions only

3. VERIFY TEST DATA
   - Check starter_code/test_data/meshes/ contains 4 OBJ files (cube, cylinder, sphere, torus)

4. BUILD PART 1 (C++)
   cd starter_code/part1_cpp
   mkdir build && cd build
   cmake ..
   make
   ./test_unwrap

   Expected: 10 tests (will fail until you implement)

5. SETUP PART 2 (Python)
   cd starter_code/part2_python
   pip install -r requirements.txt

6. PART 3 (Blender)
   - No starter code provided
   - Create your add-on from scratch
   - See ASSIGNMENT.md for requirements

================================================================================
REPOSITORY STRUCTURE
================================================================================

mixar-uv-unwrap-assignment/
├── README.txt                     ← YOU ARE HERE
├── README.md                      ← README FOR YOUR CODE
├── ASSIGNMENT.md                  ← COMPLETE REQUIREMENTS
│
└── starter_code/                  ← YOUR IMPLEMENTATION GOES HERE
    ├── part1_cpp/                 (C++ unwrapping engine)
    ├── part2_python/              (Python batch processor)
    └── meshes/                    (4 test OBJ files)

NOTE: Part 3 (Blender add-on) has no starter code - you create from scratch

================================================================================
WHAT TO IMPLEMENT
================================================================================

PART 1 - C++ ENGINE (70 points, 8-10 hours)
--------------------------------------------
Files you'll implement:
  ✏️  src/topology.cpp          - Build edge connectivity
  ✏️  src/seam_detection.cpp    - Spanning tree + angular defect
  ✏️  src/lscm.cpp              - LSCM parameterization
  ✏️  src/packing.cpp           - Island packing
  ✏️  src/unwrap.cpp            - Main orchestrator

Files provided (DO NOT MODIFY):
  ✅  include/mesh.h             - Mesh data structure
  ✅  include/topology.h         - Topology API
  ✅  src/mesh_io.cpp            - OBJ file I/O
  ✅  src/math_utils.cpp         - Vector math
  ✅  tests/test_unwrap.cpp      - Test suite

Reference materials:
  📚  reference/algorithms.md              - Algorithm descriptions
  📚  reference/lscm_math.md               - Mathematical background
  📚  reference/lscm_matrix_example.cpp    - LSCM example

PART 2 - PYTHON PROCESSOR (35 points, 4-5 hours)
-------------------------------------------------
Files you'll implement:
  ✏️  uvwrap/bindings.py        - C++ library wrapper
  ✏️  uvwrap/processor.py       - Multi-threaded batch processing
  ✏️  uvwrap/metrics.py         - Quality metrics
  ✏️  uvwrap/optimizer.py       - Parameter optimization
  ✏️  cli.py                    - Command-line interface

Reference materials:
  📚  reference/metrics_spec.md        - Exact metric formulas with Python examples

PART 3 - BLENDER ADD-ON (35 points, 6-8 hours)
-----------------------------------------------
⚠️  NO STARTER CODE PROVIDED - Build from scratch

Part 3 tests your ability to create production tools independently.
You will create a complete Blender add-on that:
  ✏️  Integrates with Blender 4.2+
  ✏️  Calls your C++ unwrapping engine via Python bindings
  ✏️  Provides UI for parameter control
  ✏️  Implements caching system
  ✏️  Supports batch processing
  ✏️  Includes seam editing tools

See ASSIGNMENT.md lines 142-193 for complete requirements.

You will create files like:
  - __init__.py (add-on registration with bl_info)
  - operators.py (unwrap, batch, seam operators)
  - panels.py (UI panel in 3D viewport)
  - core/cache.py (caching system)
  - Any other files you need for your architecture

================================================================================
DEPENDENCIES
================================================================================

PART 1 - C++:
  - CMake 3.15+
  - C++14 compiler (GCC 7+, Clang 6+, MSVC 2017+)
  - Eigen 3.3+ linear algebra library

EIGEN INSTALLATION:
  Ubuntu/Debian:  sudo apt-get install libeigen3-dev
  macOS:          brew install eigen
  Windows:        vcpkg install eigen3
  From source:    Download from eigen.tuxfamily.org

PART 2 - Python:
  - Python 3.8+
  - numpy
  - See requirements.txt for full list

PART 3 - Blender:
  - Blender 4.2+ (download from blender.org)
  - Uses Blender's bundled Python

================================================================================
DEVELOPMENT WORKFLOW
================================================================================

RECOMMENDED ORDER:
1. Part 1: C++ Engine (foundation for everything)
   - Start with topology.cpp (simplest)
   - Then seam_detection.cpp
   - Then lscm.cpp (most complex - budget 4-5 hours)
   - Finally packing.cpp

2. Part 2: Python Processor (builds on Part 1)
   - Start with bindings.py to wrap C++ library
   - Implement metrics.py
   - Add processor.py for multi-threading
   - Build cli.py
   - Add optimizer.py

3. Part 3: Blender Add-on (integrates everything)
   - Start with basic unwrap operator
   - Add UI panel
   - Implement caching
   - Add seam editing
   - Add batch processing
   - Add live preview

TESTING:
- Part 1: Run ./test_unwrap after each implementation
- Part 2: Run pytest tests/ after each module
- Part 3: Test in Blender after each feature

VALIDATION:
Run the test suite to verify your implementation:
cd starter_code/part1_cpp/build
./test_unwrap

MEMORY LEAK CHECKING:
Check your implementation for memory leaks using valgrind:

# Install valgrind (if not already installed)
# Ubuntu/Debian:
sudo apt-get install valgrind

# macOS (Homebrew):
brew install valgrind

# Run tests with valgrind
cd starter_code/part1_cpp/build
valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         ./test_unwrap

Expected output for correct implementation:
  LEAK SUMMARY:
    definitely lost: 0 bytes in 0 blocks
    indirectly lost: 0 bytes in 0 blocks

Common issues:
  "still reachable: XXX bytes" - Usually OK (global allocations)
  "definitely lost: XXX bytes" - Memory leak! Fix required.
  "invalid read/write" - Buffer overflow! Check array bounds.

================================================================================
SUBMISSION CHECKLIST
================================================================================

Before submitting, ensure you have:

PART 1:
  ☐  All tests pass (10/10)
  ☐  ALGORITHM.md explains your approach
  ☐  TEST_RESULTS.txt contains test output
  ☐  Code compiles without warnings
  ☐  No memory leaks (run valgrind)

PART 2:
  ☐  All CLI commands work
  ☐  Multi-threading functional
  ☐  Quality metrics match spec
  ☐  README.md has usage examples

PART 3:
  ☐  Add-on installs in Blender
  ☐  All operators functional
  ☐  Caching works correctly
  ☐  demo.blend included
  ☐  Screenshots included

OVERALL:
  ☐  README.md in root with build instructions
  ☐  TIME_LOG.md with time breakdown
  ☐  DESIGN_DECISIONS.md with rationale
  ☐  All code documented
  ☐  Validated against reference outputs

================================================================================
GETTING HELP
================================================================================

QUESTIONS:
- Email: ajay@mixar.app
- Subject: "Assignment Question - [Your Name]"
- Ask within first 72 hours for best response time
- We respond within 1 business day (weekdays)

DEBUGGING:
1. Check reference/ directories for examples
2. Run test suite to see specific failure messages
3. Check test output for detailed error information
4. Use printf debugging in skeleton code

COMMON ISSUES:
- "CMake can't find Eigen" → Install system package (see EIGEN INSTALLATION above)
- "Tests segfault" → Check memory management, see MEMORY LEAK CHECKING section
- "LSCM doesn't converge" → Check boundary conditions, matrix assembly
- "Blender can't import addon" → Check __init__.py has bl_info

================================================================================
RESOURCES
================================================================================

ALGORITHMS:
- starter_code/part1_cpp/reference/algorithms.md
- starter_code/part1_cpp/reference/lscm_math.md

EXAMPLES:
- All reference/ directories contain working examples
- Study these before implementing

TEST DATA:
- test_data/meshes/ - 4 test meshes for validation

================================================================================
GOOD LUCK!
================================================================================

This assignment reflects real work at Mixar. We're excited to see your solution!

Tips for success:
  ✓ Start early, ask questions
  ✓ Test incrementally
  ✓ Document as you go
  ✓ Use reference implementations

— Mixar Team

================================================================================
