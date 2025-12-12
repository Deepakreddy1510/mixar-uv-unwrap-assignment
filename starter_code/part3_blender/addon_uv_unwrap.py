bl_info = {
    "name": "MixAR UV Unwrap (Part3)",
    "author": "Auto-generated helper",
    "version": (1, 0),
    "blender": (5, 0, 0),
    "location": "View3D > Sidebar > MixAR UV",
    "description": "Unwrap using C++ binding (if present) + caching, seam editing, batch, live preview",
    "category": "UV",
}

import bpy, bmesh, time, json, os, hashlib
from bpy.props import FloatProperty, IntProperty, BoolProperty, StringProperty
from mathutils import Vector

# ---- Configuration / cache utilities ----
CACHE_DIR = os.path.join(os.path.dirname(__file__), ".part3_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def mesh_hash(obj, params: dict):
    """Hash vertex positions and faces and params to make cache key."""
    me = obj.data
    verts = [v.co.to_tuple(6) for v in me.vertices]
    faces = [tuple(p.vertices) for p in me.polygons]
    h = hashlib.sha1()
    # encode geometry
    for v in verts:
        h.update((",").join(map(lambda x: f"{x:.6f}", v)).encode())
    for f in faces:
        h.update(("_".join(map(str,f))).encode())
    # encode params
    for k in sorted(params.keys()):
        h.update(f"{k}:{params[k]}".encode())
    return h.hexdigest()

def cache_save(key: str, data: dict):
    p = os.path.join(CACHE_DIR, f"{key}.json")
    with open(p, "w") as f:
        json.dump({"ts": time.time(), "data": data}, f)

def cache_load(key: str):
    p = os.path.join(CACHE_DIR, f"{key}.json")
    if not os.path.exists(p):
        return None
    try:
        with open(p, "r") as f:
            obj = json.load(f)
            return obj.get("data")
    except:
        return None

def cache_cleanup(max_age_seconds=7*24*3600):
    # remove cache files older than max_age_seconds
    now = time.time()
    for name in os.listdir(CACHE_DIR):
        if not name.endswith(".json"):
            continue
        p = os.path.join(CACHE_DIR, name)
        try:
            ts = os.path.getmtime(p)
            if now - ts > max_age_seconds:
                os.remove(p)
        except:
            pass

# ---- binding helper: try to import your C++ binding from part2_python.bindings ----
def call_cpp_unwrap(vertices, triangles, params):
    """
    Attempt to call the C++ binding. Expected signature in bindings:
      def unwrap(vertices: List[Tuple3], triangles: List[Tuple3], params: dict) -> List[Tuple2] (per-vertex UVs)
    If binding missing or raises, return None.
    """
    try:
        # import dynamic. Adjust path if your binding module name differs.
        import importlib, sys
        # ensure repo path is in sys.path
        repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
        if repo_root not in sys.path:
            sys.path.insert(0, repo_root)
        bindings = importlib.import_module("starter_code.part2_python.bindings")
        if hasattr(bindings, "unwrap"):
            uvs = bindings.unwrap(vertices, triangles, params)
            return uvs
        elif hasattr(bindings, "unwrap_mesh"):
            uvs = bindings.unwrap_mesh(vertices, triangles, params)
            return uvs
        else:
            return None
    except Exception as e:
        print("Binding unwrap failed:", e)
        return None

# ---- mesh extraction / apply helpers ----
def extract_mesh_data(obj):
    """Return list of verts (x,y,z) and triangles (i0,i1,i2). Works on mesh objects."""
    me = obj.data
    # ensure up-to-date
    me.calc_loop_triangles()
    vertices = [tuple(v.co) for v in me.vertices]
    triangles = []
    for tri in me.loop_triangles:
        triangles.append(tuple([tri.vertices[0], tri.vertices[1], tri.vertices[2]]))
    return vertices, triangles

def apply_uvs_to_mesh(obj, uvs, per_vertex=True):
    """Apply uvs: if per_vertex True len(uvs)==len(vertices), else uv per-loop tri-based list."""
    me = obj.data
    # ensure UV layer
    if me.uv_layers.active is None:
        me.uv_layers.new(name="MixAR_UV")
    uv_layer = me.uv_layers.active.data
    me.calc_loop_triangles()

    if per_vertex:
        # map per-loop
        # uv per vertex -> set each loop's uv from vertex index
        for loop in me.loops:
            vi = loop.vertex_index
            uv_layer[loop.index].uv = Vector(uvs[vi])
    else:
        # uv per-face-vertex: assume uvs length == len(loop_triangles)*3
        idx = 0
        for tri in me.loop_triangles:
            loops = tri.loops
            for li in loops:
                uv_layer[li].uv = Vector(uvs[idx])
                idx += 1

# ---- Operators ----

class OBJECT_OT_mixar_unwrap(bpy.types.Operator):
    """Unwrap mesh using C++ binding or Blender fallback"""
    bl_idname = "object.mixar_unwrap"
    bl_label = "MixAR: Unwrap Mesh"
    bl_options = {'REGISTER', 'UNDO'}

    angle_threshold: FloatProperty(name="Angle Threshold", default=30.0)
    min_island: IntProperty(name="Min Island Size", default=1, min=1)
    margin: FloatProperty(name="Island Margin", default=0.02)

    def execute(self, context):
        obj = context.object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "Select a mesh object")
            return {'CANCELLED'}

        params = {"angle_threshold": self.angle_threshold,
                  "min_island": self.min_island,
                  "margin": self.margin}
        key = mesh_hash(obj, params)
        cached = cache_load(key)
        if cached:
            # apply cached UVs
            apply_uvs_to_mesh(obj, cached.get("uvs"), per_vertex=cached.get("per_vertex", True))
            self.report({'INFO'}, "Applied UVs from cache")
            return {'FINISHED'}

        verts, tris = extract_mesh_data(obj)
        # try C++ binding
        uvs = call_cpp_unwrap(verts, tris, params)
        if uvs:
            # assume per-vertex uvs if length matches
            per_v = (len(uvs) == len(verts))
            apply_uvs_to_mesh(obj, uvs, per_vertex=per_v)
            cache_save(key, {"uvs": uvs, "per_vertex": per_v})
            self.report({'INFO'}, "Unwrapped via C++ binding and cached")
            return {'FINISHED'}
        else:
            # fallback - use Blender Smart UV Project on a copy and copy UVs back
            bpy.ops.object.mode_set(mode='OBJECT')
            # duplicate mesh data in-memory by temporarily using Smart UV on the object
            # store selected state
            prev_mode = context.mode
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.select_all(action='SELECT')
            bpy.ops.uv.smart_project(angle_limit=self.angle_threshold, island_margin=self.margin)
            bpy.ops.object.mode_set(mode='OBJECT')
            self.report({'WARNING'}, "C++ binding not available; used Blender Smart UV Project")
            # extract resulting uvs and cache them
            me = obj.data
            me.calc_loop_triangles()
            # get per-loop uvs
            uvs_loop = [tuple(uv.uv) for uv in me.uv_layers.active.data]
            # fallback store as per-vertex by averaging per-vertex loops
            v_count = len(me.vertices)
            accum = [Vector((0.0,0.0)) for _ in range(v_count)]
            cnt = [0]*v_count
            for loop in me.loops:
                vi = loop.vertex_index
                accum[vi] += Vector(me.uv_layers.active.data[loop.index].uv)
                cnt[vi] += 1
            per_vertex_uvs = [ (accum[i]/max(1,cnt[i])) for i in range(v_count) ]
            per_vertex_list = [(uv.x, uv.y) for uv in per_vertex_uvs]
            cache_save(key, {"uvs": per_vertex_list, "per_vertex": True})
            return {'FINISHED'}

class OBJECT_OT_mixar_mark_seams(bpy.types.Operator):
    """Toggle seam on selected edges and invalidate cache"""
    bl_idname = "mesh.mixar_toggle_seam"
    bl_label = "MixAR: Toggle Seams"
    bl_options = {'REGISTER', 'UNDO'}

    make_seam: BoolProperty(name="Mark Seam", default=True)

    def execute(self, context):
        obj = context.object
        if obj is None or obj.type != 'MESH':
            self.report({'ERROR'}, "Select a mesh in Object Mode")
            return {'CANCELLED'}
        bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(obj.data)
        edges = [e for e in bm.edges if e.select]
        if not edges:
            self.report({'WARNING'}, "No edges selected")
            bpy.ops.object.mode_set(mode='OBJECT')
            return {'CANCELLED'}
        for e in edges:
            e.seam = self.make_seam
        bmesh.update_edit_mesh(obj.data)
        bpy.ops.object.mode_set(mode='OBJECT')
        # Invalidate cache entries that involve this mesh: simple strategy - clear whole cache
        for fn in os.listdir(CACHE_DIR):
            if fn.endswith(".json"):
                try:
                    os.remove(os.path.join(CACHE_DIR, fn))
                except:
                    pass
        self.report({'INFO'}, "Toggled seams and cleared cache")
        return {'FINISHED'}

class OBJECT_OT_mixar_batch_unwrap(bpy.types.Operator):
    """Batch unwrap all mesh objects in the scene"""
    bl_idname = "object.mixar_batch_unwrap"
    bl_label = "MixAR: Batch Unwrap Scene"
    bl_options = {'REGISTER'}

    angle_threshold: FloatProperty(name="Angle Threshold", default=30.0)
    min_island: IntProperty(name="Min Island Size", default=1, min=1)
    margin: FloatProperty(name="Island Margin", default=0.02)

    def execute(self, context):
        objs = [o for o in context.scene.objects if o.type=='MESH']
        total = len(objs)
        i = 0
        for obj in objs:
            i += 1
            context.view_layer.objects.active = obj
            bpy.ops.object.mixar_unwrap('INVOKE_DEFAULT',
                                       angle_threshold=self.angle_threshold,
                                       min_island=self.min_island,
                                       margin=self.margin)
        self.report({'INFO'}, f"Batch unwrapped {total} meshes")
        return {'FINISHED'}

# Modal live-preview operator with simple debouncing and FPS
class VIEW3D_OT_mixar_live_preview(bpy.types.Operator):
    """Live preview unwrap (modal)"""
    bl_idname = "view3d.mixar_live_preview"
    bl_label = "MixAR: Live Preview"
    _timer = None
    _last_run = 0.0
    _frame_times = []

    debounce_ms = 200  # 200ms

    def modal(self, context, event):
        if event.type == 'TIMER':
            now = time.time()
            dt = now - self._last_run
            if dt*1000.0 >= self.debounce_ms:
                obj = context.object
                if obj and obj.type=='MESH':
                    # call unwrap operator non-blocking
                    bpy.ops.object.mixar_unwrap('INVOKE_DEFAULT')
                self._last_run = now
                # fps tracking
                self._frame_times.append(now)
                # keep last 20 times
                if len(self._frame_times) > 20:
                    self._frame_times.pop(0)
            # draw FPS in the header by setting a window property (simple)
            if self._frame_times and len(self._frame_times) > 1:
                fps = len(self._frame_times) / (self._frame_times[-1] - self._frame_times[0] + 1e-6)
                context.window_manager.clipboard = f"MixAR Live FPS: {fps:.1f}"
        if event.type in {'ESC'}:
            self.cancel(context)
            return {'CANCELLED'}
        return {'PASS_THROUGH'}

    def invoke(self, context, event):
        wm = context.window_manager
        self._timer = wm.event_timer_add(0.1, window=context.window)  # 10 Hz
        wm.modal_handler_add(self)
        self._last_run = 0.0
        self._frame_times = []
        self.report({'INFO'}, "Live preview started. Press ESC to stop.")
        return {'RUNNING_MODAL'}

    def cancel(self, context):
        wm = context.window_manager
        if self._timer is not None:
            wm.event_timer_remove(self._timer)
        self.report({'INFO'}, "Live preview stopped.")

# ---- UI Panel ----
class VIEW3D_PT_mixar_panel(bpy.types.Panel):
    bl_label = "MixAR UV"
    bl_category = "MixAR UV"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_context = "objectmode"

    def draw(self, context):
        layout = self.layout
        col = layout.column()
        col.label(text="Unwrap / Cache")
        row = col.row()
        row.operator("object.mixar_unwrap", text="Unwrap Selected")
        row = col.row()
        row.operator("object.mixar_batch_unwrap", text="Batch Unwrap")
        row = col.row()
        row.operator("mesh.mixar_toggle_seam", text="Toggle Seam (Edit Mode)")
        row = col.row()
        row.operator("view3d.mixar_live_preview", text="Start Live Preview")
        col.separator()
        col.label(text="Cache")
        row = col.row()
        cache_files = [f for f in os.listdir(CACHE_DIR) if f.endswith(".json")]
        row.label(text=f"{len(cache_files)} entries")
        row = col.row()
        row.operator("wm.mixar_clear_cache", text="Clear Cache")

class WM_OT_mixar_clear_cache(bpy.types.Operator):
    bl_idname = "wm.mixar_clear_cache"
    bl_label = "Clear MixAR Cache"
    def execute(self, context):
        for fn in os.listdir(CACHE_DIR):
            if fn.endswith(".json"):
                try:
                    os.remove(os.path.join(CACHE_DIR, fn))
                except:
                    pass
        self.report({'INFO'}, "Cache cleared")
        return {'FINISHED'}

# ---- register ----
classes = (
    OBJECT_OT_mixar_unwrap,
    OBJECT_OT_mixar_mark_seams,
    OBJECT_OT_mixar_batch_unwrap,
    VIEW3D_OT_mixar_live_preview,
    VIEW3D_PT_mixar_panel,
    WM_OT_mixar_clear_cache,
)

def register():
    for c in classes:
        bpy.utils.register_class(c)

def unregister():
    for c in reversed(classes):
        bpy.utils.unregister_class(c)

if __name__ == "__main__":
    register()
