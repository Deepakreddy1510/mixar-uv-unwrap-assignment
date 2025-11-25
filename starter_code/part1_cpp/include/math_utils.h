/**
 * @file math_utils.h
 * @brief Math utilities for vector operations
 *
 * PROVIDED - Implementation in math_utils.cpp
 */

#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include "mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief 3D vector
 */
typedef struct {
    float x, y, z;
} Vec3;

/**
 * @brief 2D vector
 */
typedef struct {
    float x, y;
} Vec2;

/* Vector operations */
Vec3 vec3_add(Vec3 a, Vec3 b);
Vec3 vec3_sub(Vec3 a, Vec3 b);
Vec3 vec3_scale(Vec3 v, float s);
float vec3_dot(Vec3 a, Vec3 b);
Vec3 vec3_cross(Vec3 a, Vec3 b);
float vec3_length(Vec3 v);
Vec3 vec3_normalize(Vec3 v);

Vec2 vec2_add(Vec2 a, Vec2 b);
Vec2 vec2_sub(Vec2 a, Vec2 b);
float vec2_dot(Vec2 a, Vec2 b);
float vec2_length(Vec2 v);

/* Utility functions */
float clamp_float(float v, float min_val, float max_val);
float min_float(float a, float b);
float max_float(float a, float b);

/**
 * @brief Get vertex position as Vec3
 */
Vec3 get_vertex_position(const Mesh* mesh, int vertex_idx);

/**
 * @brief Compute angle at vertex in triangle
 * @param mesh Input mesh
 * @param tri_idx Triangle index
 * @param vert_idx Vertex index
 * @return Angle in radians
 */
float compute_vertex_angle_in_triangle(const Mesh* mesh,
                                       int tri_idx,
                                       int vert_idx);

#ifdef __cplusplus
}
#endif

#endif /* MATH_UTILS_H */
