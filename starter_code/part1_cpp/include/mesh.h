/**
 * @file mesh.h
 * @brief Mesh data structure
 *
 * PROVIDED - DO NOT MODIFY
 * This file defines the core mesh data structure used throughout the system.
 */

#ifndef MESH_H
#define MESH_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Triangle mesh structure with optional UVs
 *
 * Memory layout:
 * - vertices: flat array [x,y,z, x,y,z, ...]
 * - triangles: flat array [v0,v1,v2, v0,v1,v2, ...]
 * - uvs: flat array [u,v, u,v, ...] (can be NULL)
 */
typedef struct {
    float* vertices;        /**< Vertex positions (3 * num_vertices) */
    int num_vertices;       /**< Number of vertices */

    int* triangles;         /**< Triangle indices (3 * num_triangles) */
    int num_triangles;      /**< Number of triangles */

    float* uvs;             /**< UV coordinates (2 * num_vertices), can be NULL */
} Mesh;

/**
 * @brief Load mesh from OBJ file
 * @param filename Path to OBJ file
 * @return Newly allocated mesh, or NULL on error
 * @note Caller must free with free_mesh()
 */
Mesh* load_obj(const char* filename);

/**
 * @brief Save mesh to OBJ file
 * @param mesh Mesh to save
 * @param filename Output path
 * @return 0 on success, -1 on error
 */
int save_obj(const Mesh* mesh, const char* filename);

/**
 * @brief Free mesh memory
 * @param mesh Mesh to free
 */
void free_mesh(Mesh* mesh);

/**
 * @brief Allocate new mesh with same topology as input
 * @param input Source mesh
 * @return Newly allocated mesh (UVs set to NULL)
 */
Mesh* allocate_mesh_copy(const Mesh* input);

#ifdef __cplusplus
}
#endif

#endif /* MESH_H */
