/**
 * @file mesh_io.cpp
 * @brief OBJ file I/O for meshes
 *
 * PROVIDED - Complete implementation
 * Handles loading and saving OBJ files with UVs
 */

#include "mesh.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>

Mesh* load_obj(const char* filename) {
    FILE* f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Cannot open file: %s\n", filename);
        return NULL;
    }

    std::vector<float> vertices;
    std::vector<int> triangles;
    std::vector<float> uvs_temp;
    bool has_uvs = false;

    char line[256];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == 'v' && line[1] == ' ') {
            // Vertex
            float x, y, z;
            if (sscanf(line, "v %f %f %f", &x, &y, &z) == 3) {
                vertices.push_back(x);
                vertices.push_back(y);
                vertices.push_back(z);
            }
        } else if (line[0] == 'v' && line[1] == 't') {
            // UV coordinate
            float u, v;
            if (sscanf(line, "vt %f %f", &u, &v) == 2) {
                uvs_temp.push_back(u);
                uvs_temp.push_back(v);
                has_uvs = true;
            }
        } else if (line[0] == 'f' && line[1] == ' ') {
            // Face - supports multiple formats:
            // - Triangles: f v1 v2 v3
            // - Triangles with UVs: f v1/vt1 v2/vt2 v3/vt3
            // - Triangles with UVs and normals: f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
            // - Quads: f v1 v2 v3 v4 (converted to two triangles)
            // - Quads with UVs/normals: f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3 v4/vt4/vn4

            int v[4], vt[4], vn[4];
            int num_vertices = vertices.size() / 3;
            int num_parsed = 0;

            // Try v/vt/vn format (most common in Blender)
            if (sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",
                      &v[0], &vt[0], &vn[0], &v[1], &vt[1], &vn[1],
                      &v[2], &vt[2], &vn[2], &v[3], &vt[3], &vn[3]) == 12) {
                num_parsed = 4;  // Quad with UVs and normals
            } else if (sscanf(line, "f %d/%d/%d %d/%d/%d %d/%d/%d",
                             &v[0], &vt[0], &vn[0], &v[1], &vt[1], &vn[1],
                             &v[2], &vt[2], &vn[2]) == 9) {
                num_parsed = 3;  // Triangle with UVs and normals
            } else if (sscanf(line, "f %d/%d %d/%d %d/%d %d/%d",
                             &v[0], &vt[0], &v[1], &vt[1], &v[2], &vt[2], &v[3], &vt[3]) == 8) {
                num_parsed = 4;  // Quad with UVs
            } else if (sscanf(line, "f %d/%d %d/%d %d/%d",
                             &v[0], &vt[0], &v[1], &vt[1], &v[2], &vt[2]) == 6) {
                num_parsed = 3;  // Triangle with UVs
            } else if (sscanf(line, "f %d %d %d %d", &v[0], &v[1], &v[2], &v[3]) == 4) {
                num_parsed = 4;  // Quad without UVs
            } else if (sscanf(line, "f %d %d %d", &v[0], &v[1], &v[2]) == 3) {
                num_parsed = 3;  // Triangle without UVs
            }

            if (num_parsed >= 3) {
                // Validate vertex indices (OBJ is 1-indexed)
                bool valid = true;
                for (int i = 0; i < num_parsed; i++) {
                    if (v[i] < 1 || v[i] > num_vertices) {
                        fprintf(stderr, "Error: Invalid vertex index %d in face (valid range: 1-%d)\n",
                               v[i], num_vertices);
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    // Add first triangle (v0, v1, v2)
                    triangles.push_back(v[0] - 1);
                    triangles.push_back(v[1] - 1);
                    triangles.push_back(v[2] - 1);

                    // If quad, add second triangle (v0, v2, v3)
                    if (num_parsed == 4) {
                        triangles.push_back(v[0] - 1);
                        triangles.push_back(v[2] - 1);
                        triangles.push_back(v[3] - 1);
                    }
                }
            }
        }
    }

    fclose(f);

    if (vertices.empty() || triangles.empty()) {
        fprintf(stderr, "Failed to parse OBJ file: %s\n", filename);
        return NULL;
    }

    Mesh* mesh = (Mesh*)malloc(sizeof(Mesh));

    mesh->num_vertices = vertices.size() / 3;
    mesh->vertices = (float*)malloc(vertices.size() * sizeof(float));
    memcpy(mesh->vertices, vertices.data(), vertices.size() * sizeof(float));

    mesh->num_triangles = triangles.size() / 3;
    mesh->triangles = (int*)malloc(triangles.size() * sizeof(int));
    memcpy(mesh->triangles, triangles.data(), triangles.size() * sizeof(int));

    if (has_uvs) {
        size_t expected_uv_count = (vertices.size() / 3) * 2;
        if (uvs_temp.size() == expected_uv_count) {
            mesh->uvs = (float*)malloc(uvs_temp.size() * sizeof(float));
            memcpy(mesh->uvs, uvs_temp.data(), uvs_temp.size() * sizeof(float));
        } else {
            fprintf(stderr,
                    "Warning: UV count mismatch in %s\n"
                    "  Expected: %zu UVs (%zu vertices)\n"
                    "  Found:    %zu UVs\n"
                    "  UVs will be ignored.\n",
                    filename, expected_uv_count / 2, vertices.size() / 3,
                    uvs_temp.size() / 2);
            mesh->uvs = NULL;
        }
    } else {
        mesh->uvs = NULL;
    }

    printf("Loaded %s: %d vertices, %d triangles\n",
           filename, mesh->num_vertices, mesh->num_triangles);

    return mesh;
}

int save_obj(const Mesh* mesh, const char* filename) {
    if (!mesh) return -1;

    FILE* f = fopen(filename, "w");
    if (!f) {
        fprintf(stderr, "Cannot write file: %s\n", filename);
        return -1;
    }

    // Write vertices
    for (int i = 0; i < mesh->num_vertices; i++) {
        fprintf(f, "v %f %f %f\n",
                mesh->vertices[i*3],
                mesh->vertices[i*3+1],
                mesh->vertices[i*3+2]);
    }

    // Write UVs if present
    if (mesh->uvs) {
        for (int i = 0; i < mesh->num_vertices; i++) {
            fprintf(f, "vt %f %f\n",
                    mesh->uvs[i*2],
                    mesh->uvs[i*2+1]);
        }
    }

    // Write faces
    for (int i = 0; i < mesh->num_triangles; i++) {
        int v0 = mesh->triangles[i*3] + 1;
        int v1 = mesh->triangles[i*3+1] + 1;
        int v2 = mesh->triangles[i*3+2] + 1;

        if (mesh->uvs) {
            fprintf(f, "f %d/%d %d/%d %d/%d\n",
                    v0, v0, v1, v1, v2, v2);
        } else {
            fprintf(f, "f %d %d %d\n", v0, v1, v2);
        }
    }

    fclose(f);
    printf("Saved %s\n", filename);
    return 0;
}

void free_mesh(Mesh* mesh) {
    if (!mesh) return;

    if (mesh->vertices) free(mesh->vertices);
    if (mesh->triangles) free(mesh->triangles);
    if (mesh->uvs) free(mesh->uvs);
    free(mesh);
}

Mesh* allocate_mesh_copy(const Mesh* input) {
    if (!input) return NULL;

    Mesh* mesh = (Mesh*)malloc(sizeof(Mesh));

    mesh->num_vertices = input->num_vertices;
    mesh->vertices = (float*)malloc(input->num_vertices * 3 * sizeof(float));
    memcpy(mesh->vertices, input->vertices, input->num_vertices * 3 * sizeof(float));

    mesh->num_triangles = input->num_triangles;
    mesh->triangles = (int*)malloc(input->num_triangles * 3 * sizeof(int));
    memcpy(mesh->triangles, input->triangles, input->num_triangles * 3 * sizeof(int));

    mesh->uvs = NULL;

    return mesh;
}
