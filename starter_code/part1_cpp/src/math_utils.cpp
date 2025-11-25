/**
 * @file math_utils.cpp
 * @brief Math utilities implementation
 *
 * PROVIDED - Complete implementation of vector math
 */

#include "math_utils.h"
#include "mesh.h"
#include <math.h>

/* Vector3 operations */
Vec3 vec3_add(Vec3 a, Vec3 b) {
    Vec3 result = {a.x + b.x, a.y + b.y, a.z + b.z};
    return result;
}

Vec3 vec3_sub(Vec3 a, Vec3 b) {
    Vec3 result = {a.x - b.x, a.y - b.y, a.z - b.z};
    return result;
}

Vec3 vec3_scale(Vec3 v, float s) {
    Vec3 result = {v.x * s, v.y * s, v.z * s};
    return result;
}

float vec3_dot(Vec3 a, Vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 vec3_cross(Vec3 a, Vec3 b) {
    Vec3 result;
    result.x = a.y * b.z - a.z * b.y;
    result.y = a.z * b.x - a.x * b.z;
    result.z = a.x * b.y - a.y * b.x;
    return result;
}

float vec3_length(Vec3 v) {
    return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec3 vec3_normalize(Vec3 v) {
    float len = vec3_length(v);
    if (len < 1e-8f) {
        Vec3 zero = {0, 0, 0};
        return zero;
    }
    return vec3_scale(v, 1.0f / len);
}

/* Vector2 operations */
Vec2 vec2_add(Vec2 a, Vec2 b) {
    Vec2 result = {a.x + b.x, a.y + b.y};
    return result;
}

Vec2 vec2_sub(Vec2 a, Vec2 b) {
    Vec2 result = {a.x - b.x, a.y - b.y};
    return result;
}

float vec2_dot(Vec2 a, Vec2 b) {
    return a.x * b.x + a.y * b.y;
}

float vec2_length(Vec2 v) {
    return sqrtf(v.x * v.x + v.y * v.y);
}

/* Utility functions */
float clamp_float(float v, float min_val, float max_val) {
    if (v < min_val) return min_val;
    if (v > max_val) return max_val;
    return v;
}

float min_float(float a, float b) {
    return (a < b) ? a : b;
}

float max_float(float a, float b) {
    return (a > b) ? a : b;
}

Vec3 get_vertex_position(const Mesh* mesh, int vertex_idx) {
    Vec3 pos;
    pos.x = mesh->vertices[vertex_idx * 3];
    pos.y = mesh->vertices[vertex_idx * 3 + 1];
    pos.z = mesh->vertices[vertex_idx * 3 + 2];
    return pos;
}

float compute_vertex_angle_in_triangle(const Mesh* mesh,
                                       int tri_idx,
                                       int vert_idx) {
    int v0 = mesh->triangles[tri_idx * 3];
    int v1 = mesh->triangles[tri_idx * 3 + 1];
    int v2 = mesh->triangles[tri_idx * 3 + 2];

    Vec3 p, p1, p2;

    // Find which vertex is vert_idx and get the other two
    if (v0 == vert_idx) {
        p = get_vertex_position(mesh, v0);
        p1 = get_vertex_position(mesh, v1);
        p2 = get_vertex_position(mesh, v2);
    } else if (v1 == vert_idx) {
        p = get_vertex_position(mesh, v1);
        p1 = get_vertex_position(mesh, v2);
        p2 = get_vertex_position(mesh, v0);
    } else if (v2 == vert_idx) {
        p = get_vertex_position(mesh, v2);
        p1 = get_vertex_position(mesh, v0);
        p2 = get_vertex_position(mesh, v1);
    } else {
        return 0.0f;  // Vertex not in triangle
    }

    // Compute angle using dot product
    Vec3 e1 = vec3_normalize(vec3_sub(p1, p));
    Vec3 e2 = vec3_normalize(vec3_sub(p2, p));

    float cos_angle = vec3_dot(e1, e2);
    cos_angle = clamp_float(cos_angle, -1.0f, 1.0f);

    float angle = acosf(cos_angle);

    return angle;
}
