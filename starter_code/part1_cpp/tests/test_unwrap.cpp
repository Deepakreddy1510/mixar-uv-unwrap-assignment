/**
 * @file test_unwrap.cpp
 * @brief Test suite for UV unwrapping
 *
 * PROVIDED - Complete test suite
 * Run this to validate your implementation
 */

#include "mesh.h"
#include "topology.h"
#include "unwrap.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_DATA_DIR "../../test_data/meshes/"

int tests_passed = 0;
int tests_failed = 0;

void test_topology(const char* mesh_name, int expected_v, int expected_e, int expected_f) {
    printf("[TEST] Topology - %s...", mesh_name);

    char filename[256];
    snprintf(filename, sizeof(filename), "%s%s", TEST_DATA_DIR, mesh_name);

    Mesh* mesh = load_obj(filename);
    if (!mesh) {
        printf(" FAIL (could not load)\n");
        tests_failed++;
        return;
    }

    TopologyInfo* topo = build_topology(mesh);
    if (!topo) {
        printf(" FAIL (topology building failed)\n");
        tests_failed++;
        free_mesh(mesh);
        return;
    }

    int v = mesh->num_vertices;
    int e = topo->num_edges;
    int f = mesh->num_triangles;

    if (v != expected_v || e != expected_e || f != expected_f) {
        printf(" FAIL\n");
        printf("  Expected: V=%d, E=%d, F=%d\n", expected_v, expected_e, expected_f);
        printf("  Got:      V=%d, E=%d, F=%d\n", v, e, f);
        tests_failed++;
    } else {
        int euler = v - e + f;
        printf(" PASS (Euler: %d)\n", euler);
        tests_passed++;
    }

    free_topology(topo);
    free_mesh(mesh);
}

void test_seams(const char* mesh_name, int min_seams, int max_seams) {
    printf("[TEST] Seam Detection - %s...", mesh_name);

    char filename[256];
    snprintf(filename, sizeof(filename), "%s%s", TEST_DATA_DIR, mesh_name);

    Mesh* mesh = load_obj(filename);
    if (!mesh) {
        printf(" FAIL (could not load)\n");
        tests_failed++;
        return;
    }

    TopologyInfo* topo = build_topology(mesh);
    if (!topo) {
        printf(" FAIL (topology failed)\n");
        tests_failed++;
        free_mesh(mesh);
        return;
    }

    int num_seams;
    int* seams = detect_seams(mesh, topo, 30.0f, &num_seams);

    if (!seams) {
        printf(" FAIL (seam detection failed)\n");
        tests_failed++;
    } else if (num_seams < min_seams || num_seams > max_seams) {
        printf(" FAIL\n");
        printf("  Expected: %d-%d seams\n", min_seams, max_seams);
        printf("  Got:      %d seams\n", num_seams);
        tests_failed++;
    } else {
        printf(" PASS (%d seams)\n", num_seams);
        tests_passed++;
    }

    if (seams) free(seams);
    free_topology(topo);
    free_mesh(mesh);
}

void test_unwrap(const char* mesh_name, float max_stretch_threshold) {
    printf("[TEST] Unwrap - %s...", mesh_name);

    char filename[256];
    snprintf(filename, sizeof(filename), "%s%s", TEST_DATA_DIR, mesh_name);

    Mesh* mesh = load_obj(filename);
    if (!mesh) {
        printf(" FAIL (could not load)\n");
        tests_failed++;
        return;
    }

    UnwrapParams params;
    params.angle_threshold = 30.0f;
    params.min_island_faces = 5;
    params.pack_islands = 1;
    params.island_margin = 0.02f;

    UnwrapResult* result;
    Mesh* unwrapped = unwrap_mesh(mesh, &params, &result);

    if (!unwrapped || !result) {
        printf(" FAIL (unwrapping failed)\n");
        tests_failed++;
        free_mesh(mesh);
        return;
    }

    if (!unwrapped->uvs) {
        printf(" FAIL (no UVs generated)\n");
        tests_failed++;
        free_unwrap_result(result);
        free_mesh(unwrapped);
        free_mesh(mesh);
        return;
    }

    // Check quality
    float stretch = result->max_stretch;

    if (stretch > max_stretch_threshold) {
        printf(" FAIL (stretch=%.2f > %.2f)\n", stretch, max_stretch_threshold);
        tests_failed++;
    } else {
        printf(" PASS (islands=%d, stretch=%.2f)\n", result->num_islands, stretch);
        tests_passed++;
    }

    free_unwrap_result(result);
    free_mesh(unwrapped);
    free_mesh(mesh);
}

int main() {
    printf("\n");
    printf("========================================\n");
    printf("UV Unwrapping Test Suite\n");
    printf("========================================\n\n");

    // Topology tests
    test_topology("01_cube.obj", 8, 18, 12);
    test_topology("03_sphere.obj", 42, 120, 80);

    // Seam detection tests
    // Basic spanning tree should produce minimum seams
    // Angular defect refinement may add 2-4 additional seams
    test_seams("01_cube.obj", 7, 11);           // Basic: 7, refined: 7-11
    test_seams("03_sphere.obj", 1, 5);          // Sphere needs more seams due to curvature
    test_seams("02_cylinder.obj", 1, 3);        // Cylinder: 1-2 seams typically

    // Full unwrap tests
    test_unwrap("01_cube.obj", 2.0f);           // Allow up to 2.0 stretch
    test_unwrap("03_sphere.obj", 2.0f);
    test_unwrap("02_cylinder.obj", 1.5f);       // Cylinder should be better

    printf("\n");
    printf("========================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("========================================\n\n");

    return (tests_failed == 0) ? 0 : 1;
}
