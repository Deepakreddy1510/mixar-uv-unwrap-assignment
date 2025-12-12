# Part 3 — Blender Add-on (UV Unwrap)

This part implements a Blender add-on that integrates the custom UV unwrapping pipeline
developed in Part 1 (C++) and Part 2 (Python). The add-on exposes a UI panel inside Blender
and supports the following features:

### ✔ Main Unwrap Operator
- Extracts mesh data directly from Blender
- Calls the C++ unwrap library through Python bindings when available
- Falls back to Blender’s Smart UV Project
- Applies UVs directly to the active mesh

### ✔ Caching System
- Hashes mesh geometry + parameters
- Stores cached unwrap results as JSON files
- Automatically reuses cached UVs
- Supports manual “Clear Cache”

### ✔ UI Panel
Located under **View → Sidebar (N) → MixAR UV**:
- Unwrap Selected
- Batch Unwrap (unwrap all mesh objects)
- Toggle Seam (Edit Mode)
- Start Live Preview
- Cache status indicator

### ✔ Seam Editing
- Select edges in Edit Mode
- Use "Toggle Seam" to mark or clear seams
- Cache invalidates automatically

### ✔ Batch Processing
- Processes all mesh objects in the scene
- Shows operator info

### ✔ Live Preview
- Modal operator refreshes unwrap at intervals
- Debouncing (200 ms)
- FPS output through clipboard text

### Demo Files
- `demo.blend` (optional scene file)
- UV layout exported: `screenshots/uv_layout.png`
- UI screenshots: `screenshots/panel.png`

