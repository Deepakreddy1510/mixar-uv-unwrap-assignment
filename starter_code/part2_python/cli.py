"""
Command-line interface for UV unwrapping

Implemented commands:
- unwrap: Unwrap single mesh
- batch: Batch process directory (multi-threaded)
- optimize: Find optimal parameters (simple grid search)
- analyze: Analyze mesh quality (computes metrics)
"""

import argparse
import sys
import os
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import multiprocessing
from time import time

# third-party / local imports
try:
    from uvwrap import bindings
    import uvwrap.metrics as metrics
    import numpy as np
except Exception as e:
    # Defer import errors to runtime with clearer message
    bindings = None
    metrics = None
    np = None
    _IMPORT_ERROR = e
else:
    _IMPORT_ERROR = None

# Helpers ---------------------------------------------------------------------
def ensure_bindings():
    if _IMPORT_ERROR:
        raise RuntimeError(
            "Failed to import uvwrap bindings/metrics. Original error:\n" + repr(_IMPORT_ERROR)
        )

def pretty_pct(x):
    return f"{x*100:.2f}%"


def format_metrics(mdict):
    return (
        f"Stretch(max)  : {mdict.get('stretch_max', mdict.get('max_stretch', 0.0)):.6g}\n"
        f"Coverage      : {pretty_pct(mdict.get('coverage', 0.0))}\n"
        f"Angle (deg)   : {mdict.get('angle_dist_deg', mdict.get('angle_distortion_deg', 0.0)):.6g}"
    )

# Unwrap command --------------------------------------------------------------
def cmd_unwrap(args):
    """
    Unwrap single mesh and save result
    """
    ensure_bindings()

    inp = Path(args.input)
    outp = Path(args.output)

    if not inp.exists():
        print("Input file not found:", inp)
        return 1

    print(f"Unwrapping: {inp} -> {outp}")
    params = {
        "angle_threshold": float(args.angle_threshold),
        "min_island_faces": int(args.min_island),
        "pack_islands": bool(args.pack),
        "island_margin": float(args.margin),
    }

    t0 = time()
    try:
        # bindings.unwrap accepts either path or mesh object; it returns (mesh, metrics_dict)
        result_mesh, result_metrics = bindings.unwrap(str(inp), params=params)
    except TypeError:
        # fallback: maybe bindings.unwrap expects only (path) -> (mesh, metrics)
        result_mesh, result_metrics = bindings.unwrap(str(inp))
    except Exception as e:
        print("Error during unwrap:", e)
        return 1
    dt = time() - t0

    # Save
    try:
        bindings.save_mesh(result_mesh, str(outp))
    except Exception as e:
        print("Warning: failed to save output OBJ:", e)

    # Compute / normalize metrics representation
    # Some variants return different metric field names; normalize them.
    m = {}
    if isinstance(result_metrics, dict):
        m["stretch_max"] = float(result_metrics.get("max_stretch", result_metrics.get("stretch_max", 0.0)))
        m["coverage"] = float(result_metrics.get("coverage", 0.0))
        m["angle_dist_rad"] = float(result_metrics.get("angle_distortion", result_metrics.get("angle_dist_rad", 0.0)))
        m["angle_dist_deg"] = float(np.degrees(m["angle_dist_rad"])) if np is not None else m["angle_dist_rad"]
    else:
        # fallback compute from produced mesh if metrics not returned
        try:
            uvs = result_mesh.uvs
            tris = result_mesh.triangles
            m["stretch_max"] = metrics.compute_stretch(result_mesh, uvs)
            m["coverage"] = metrics.compute_coverage(uvs, tris, resolution=256)
            m["angle_dist_rad"] = metrics.compute_angle_distortion(result_mesh, uvs)
            m["angle_dist_deg"] = float(np.degrees(m["angle_dist_rad"]))
        except Exception:
            m = {"stretch_max": 0.0, "coverage": 0.0, "angle_dist_rad": 0.0, "angle_dist_deg": 0.0}

    print("✓ Completed in {:.2f}s".format(dt))
    print(format_metrics(m))
    return 0

# Analyze command -------------------------------------------------------------
def cmd_analyze(args):
    """
    Analyze mesh quality metrics and print JSON report
    """
    ensure_bindings()

    inp = Path(args.input)
    if not inp.exists():
        print("Input file not found:", inp)
        return 1

    try:
        mesh = bindings.load_mesh(str(inp))
        if mesh is None:
            raise RuntimeError("bindings.load_mesh returned None")
    except Exception as e:
        print("Failed to load mesh via bindings:", e)
        return 1

    # If no UVs, try to unwrap with default params
    if mesh.uvs is None:
        print("Mesh has no UVs — unwrapping with default parameters...")
        mesh, _ = bindings.unwrap(str(inp))

    uvs = mesh.uvs
    tris = mesh.triangles

    # Compute metrics
    try:
        stretch = metrics.compute_stretch(mesh, uvs)
        coverage = metrics.compute_coverage(uvs, tris, resolution=512)
        angle_rad = metrics.compute_angle_distortion(mesh, uvs)
        angle_deg = float(np.degrees(angle_rad))
    except Exception as e:
        print("Error computing metrics:", e)
        return 1

    report = {
        "input": str(inp),
        "num_vertices": int(mesh.num_vertices),
        "num_triangles": int(mesh.num_triangles),
        "metrics": {
            "stretch_max": float(stretch),
            "coverage": float(coverage),
            "angle_distortion_rad": float(angle_rad),
            "angle_distortion_deg": float(angle_deg),
        }
    }

    print(json.dumps(report, indent=2))
    return 0

# Batch command ---------------------------------------------------------------
def _process_one_file(path, out_dir, params):
    """
    Worker: unwrap one file and compute metrics.
    Returns (input_path_str, success_bool, result_dict_or_error)
    """
    try:
        mesh, mm = bindings.unwrap(str(path), params=params)
    except TypeError:
        mesh, mm = bindings.unwrap(str(path))
    except Exception as e:
        return (str(path), False, {"error": str(e)})

    # build output name
    outf = out_dir / path.name
    try:
        bindings.save_mesh(mesh, str(outf))
    except Exception as e:
        # non-fatal: still compute metrics
        pass

    # normalize metrics
    if isinstance(mm, dict):
        res = {
            "stretch_max": float(mm.get("max_stretch", mm.get("stretch_max", 0.0))),
            "coverage": float(mm.get("coverage", 0.0)),
            "angle_dist_rad": float(mm.get("angle_distortion", mm.get("angle_dist_rad", 0.0))),
        }
        res["angle_dist_deg"] = float(np.degrees(res["angle_dist_rad"])) if np is not None else res["angle_dist_rad"]
    else:
        # compute from mesh
        uvs = mesh.uvs
        tris = mesh.triangles
        try:
            res = {
                "stretch_max": metrics.compute_stretch(mesh, uvs),
                "coverage": metrics.compute_coverage(uvs, tris, resolution=256),
                "angle_dist_rad": metrics.compute_angle_distortion(mesh, uvs),
            }
            res["angle_dist_deg"] = float(np.degrees(res["angle_dist_rad"]))
        except Exception as e:
            return (str(path), False, {"error": str(e)})

    return (str(path), True, res)

def cmd_batch(args):
    """
    Batch unwrap all OBJ files in input_dir and write outputs to output_dir.
    Saves optional JSON report if --report is given.
    """
    ensure_bindings()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if not input_dir.exists() or not input_dir.is_dir():
        print("Input directory not found:", input_dir)
        return 1
    output_dir.mkdir(parents=True, exist_ok=True)

    # find obj files
    obj_files = sorted([p for p in input_dir.glob("*.obj")])
    if not obj_files:
        print("No .obj files found in", input_dir)
        return 1

    # prepare params to pass to each unwrap
    params = {
        "angle_threshold": float(args.angle_threshold),
        "min_island_faces": 5,
        "pack_islands": True,
        "island_margin": 0.02,
    }

    # threads
    num_threads = int(args.threads) if args.threads else multiprocessing.cpu_count()
    print(f"Processing {len(obj_files)} files with {num_threads} threads...")

    results = {}
    with ThreadPoolExecutor(max_workers=num_threads) as ex:
        futures = {ex.submit(_process_one_file, p, output_dir, params): p for p in obj_files}
        try:
            # show progress using simple textual updates
            from tqdm import tqdm
            for fut in tqdm(as_completed(futures), total=len(futures), unit="file"):
                p = futures[fut]
                try:
                    inp, ok, res = fut.result()
                except Exception as e:
                    inp = str(p)
                    ok = False
                    res = {"error": str(e)}
                results[inp] = {"ok": ok, "result": res}
        except KeyboardInterrupt:
            print("Interrupted by user; shutting down...")
            ex.shutdown(wait=False)
            return 1

    # Print summary
    success = sum(1 for v in results.values() if v["ok"])
    failed = len(results) - success
    print(f"Done. Success: {success}; Failed: {failed}")

    if args.report:
        try:
            with open(args.report, "w") as f:
                json.dump(results, f, indent=2)
            print("Saved report to", args.report)
        except Exception as e:
            print("Failed to save report:", e)

    return 0

# Optimize command ------------------------------------------------------------
def cmd_optimize(args):
    """
    Simple grid-search optimize parameters for a single mesh.
    Minimizes the chosen metric (stretch by default).
    """
    ensure_bindings()

    inp = Path(args.input)
    if not inp.exists():
        print("Input not found:", inp)
        return 1

    # Basic parameter grid (tunable)
    angle_grid = [10.0, 20.0, 30.0, 40.0]
    min_island_grid = [1, 5, 10]
    pack_options = [True]  # keep packing enabled by default

    best = None
    best_params = None

    print(f"Running grid search ({len(angle_grid) * len(min_island_grid)} combinations) optimizing '{args.metric}'...")

    for a in angle_grid:
        for mfaces in min_island_grid:
            for pack in pack_options:
                params = {
                    "angle_threshold": float(a),
                    "min_island_faces": int(mfaces),
                    "pack_islands": bool(pack),
                    "island_margin": 0.02,
                }
                try:
                    mesh, mm = bindings.unwrap(str(inp), params=params)
                except TypeError:
                    mesh, mm = bindings.unwrap(str(inp))
                except Exception as e:
                    print("Run failed for params", params, "error:", e)
                    continue

                # normalize metric value
                metric_value = None
                if args.metric == "stretch":
                    metric_value = float(mm.get("max_stretch", mm.get("stretch_max", np.nan)))
                elif args.metric == "coverage":
                    metric_value = -float(mm.get("coverage", 0.0))  # maximize coverage => minimize negative
                elif args.metric == "angle_distortion":
                    metric_value = float(mm.get("angle_distortion", mm.get("angle_dist_rad", np.nan)))

                # handle missing values
                if metric_value is None or (isinstance(metric_value, float) and np.isnan(metric_value)):
                    continue

                if (best is None) or (metric_value < best):
                    best = metric_value
                    best_params = dict(params)

                print(f"params a={a}, min_faces={mfaces}, pack={pack} -> metric={metric_value:.6g}")

    if best_params is None:
        print("Grid search failed (no successful runs).")
        return 1

    print("Best params:", best_params)
    if args.save_params:
        try:
            with open(args.save_params, "w") as f:
                json.dump(best_params, f, indent=2)
            print("Saved best params to", args.save_params)
        except Exception as e:
            print("Failed to save params:", e)

    if args.output:
        # produce final output with best params
        try:
            mesh, mm = bindings.unwrap(str(inp), params=best_params)
            bindings.save_mesh(mesh, args.output)
            print("Saved unwrapped OBJ with best params to", args.output)
        except Exception as e:
            print("Failed to produce output with best params:", e)

    return 0

# Main ------------------------------------------------------------------------
def main():
    """
    Main CLI entry point
    """
    parser = argparse.ArgumentParser(
        description='UV Unwrapping Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Unwrap single mesh
  python cli.py unwrap input.obj output.obj --angle-threshold 30

  # Batch process
  python cli.py batch meshes/ output/ --threads 8 --report metrics.json

  # Optimize parameters
  python cli.py optimize mesh.obj --metric stretch --save-params best.json

  # Analyze quality
  python cli.py analyze mesh.obj
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Commands')

    # Unwrap command
    unwrap_parser = subparsers.add_parser('unwrap', help='Unwrap single mesh')
    unwrap_parser.add_argument('input', help='Input OBJ file')
    unwrap_parser.add_argument('output', help='Output OBJ file')
    unwrap_parser.add_argument('--angle-threshold', type=float, default=30.0,
                              help='Angle threshold in degrees (default: 30)')
    unwrap_parser.add_argument('--min-island', type=int, default=10,
                              help='Minimum island size in faces (default: 10)')
    unwrap_parser.add_argument('--no-pack', action='store_false', dest='pack',
                              default=True,
                              help='Disable island packing (default: enabled)')
    unwrap_parser.add_argument('--margin', type=float, default=0.02,
                              help='Island margin (default: 0.02)')

    # Batch command
    batch_parser = subparsers.add_parser('batch', help='Batch process directory')
    batch_parser.add_argument('input_dir', help='Input directory')
    batch_parser.add_argument('output_dir', help='Output directory')
    batch_parser.add_argument('--threads', type=int, default=None,
                             help='Number of threads (default: CPU count)')
    batch_parser.add_argument('--angle-threshold', type=float, default=30.0,
                             help='Angle threshold in degrees')
    batch_parser.add_argument('--report', help='Save metrics to JSON file')

    # Optimize command
    opt_parser = subparsers.add_parser('optimize', help='Optimize parameters')
    opt_parser.add_argument('input', help='Input OBJ file')
    opt_parser.add_argument('--output', help='Output OBJ file with best params')
    opt_parser.add_argument('--metric', choices=['stretch', 'coverage', 'angle_distortion'],
                           default='stretch', help='Metric to optimize')
    opt_parser.add_argument('--save-params', help='Save best params to JSON')

    # Analyze command
    analyze_parser = subparsers.add_parser('analyze', help='Analyze mesh quality')
    analyze_parser.add_argument('input', help='Input OBJ file with UVs')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Route to command handlers
    try:
        if args.command == 'unwrap':
            return cmd_unwrap(args)
        elif args.command == 'batch':
            return cmd_batch(args)
        elif args.command == 'optimize':
            return cmd_optimize(args)
        elif args.command == 'analyze':
            return cmd_analyze(args)
    except Exception as e:
        print("Unhandled error in CLI:", e)
        return 2

    return 0


if __name__ == '__main__':
    sys.exit(main())

