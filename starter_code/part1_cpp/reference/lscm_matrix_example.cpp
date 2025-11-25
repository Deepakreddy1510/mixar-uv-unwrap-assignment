/**
 * @file lscm_matrix_example.cpp
 * @brief Complete example of LSCM matrix assembly for ONE triangle
 *
 * REFERENCE IMPLEMENTATION
 *
 * This shows exactly how to add LSCM energy contribution for a single triangle.
 * You'll need to loop over all triangles and accumulate contributions.
 */

#include <Eigen/Sparse>
#include <cmath>
#include <iostream>

struct Vec3 {
    double x, y, z;
    Vec3(double x=0, double y=0, double z=0) : x(x), y(y), z(z) {}
    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }
};

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}

double length(const Vec3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vec3 normalize(const Vec3& v) {
    double len = length(v);
    if (len < 1e-10) return Vec3(0, 0, 0);
    return Vec3(v.x / len, v.y / len, v.z / len);
}

/**
 * Add LSCM contribution for a single triangle to the sparse matrix
 *
 * @param triplets List of matrix entries to append to
 * @param v0, v1, v2 LOCAL vertex indices in the island
 * @param p0, p1, p2 3D positions of the triangle vertices
 */
void add_triangle_lscm_contribution(
    std::vector<Eigen::Triplet<double>>& triplets,
    int v0, int v1, int v2,
    const Vec3& p0, const Vec3& p1, const Vec3& p2)
{
    // STEP 1: Project triangle to its plane
    //
    // Create local 2D coordinate system in the triangle's plane

    Vec3 e1 = p1 - p0;
    Vec3 e2 = p2 - p0;

    // Normal to triangle
    Vec3 normal = normalize(cross(e1, e2));

    // Create orthonormal basis in triangle's plane
    Vec3 u_axis = normalize(e1);
    Vec3 v_axis = cross(normal, u_axis);

    // Local 2D coordinates of triangle vertices
    double q0_x = 0.0;
    double q0_y = 0.0;

    double q1_x = dot(e1, u_axis);
    double q1_y = dot(e1, v_axis);

    double q2_x = dot(e2, u_axis);
    double q2_y = dot(e2, v_axis);

    // STEP 2: Compute triangle area
    //
    // This weights the contribution of this triangle

    double area = 0.5 * std::abs(q1_x * q2_y - q1_y * q2_x);

    if (area < 1e-10) {
        // Degenerate triangle, skip
        return;
    }

    // STEP 3: Add LSCM energy terms
    //
    // For each edge of the triangle, add contribution to the matrix
    //
    // The LSCM energy for a triangle is:
    //   E = area * ||∇u - R_90°(∇v)||²
    //
    // Where ∇u and ∇v are gradients of u and v in the triangle's plane,
    // and R_90° is a 90-degree rotation.
    //
    // This translates to matrix entries:

    // Edge v0 → v1
    {
        double dx = q1_x - q0_x;
        double dy = q1_y - q0_y;

        // Contribution to matrix
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v1,    area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v1+1,  area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v1,    area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v1+1,  area * (-dx)));

        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v0,    -area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v0,   2*v0+1,  -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v0,    -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v0+1, 2*v0+1,  -area * (-dx)));
    }

    // Edge v1 → v2
    {
        double dx = q2_x - q1_x;
        double dy = q2_y - q1_y;

        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v2,    area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v2+1,  area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v2,    area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v2+1,  area * (-dx)));

        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v1,    -area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v1,   2*v1+1,  -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v1,    -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v1+1, 2*v1+1,  -area * (-dx)));
    }

    // Edge v2 → v0
    {
        double dx = q0_x - q2_x;
        double dy = q0_y - q2_y;

        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v0,    area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v0+1,  area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v0,    area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v0+1,  area * (-dx)));

        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v2,    -area * dx));
        triplets.push_back(Eigen::Triplet<double>(2*v2,   2*v2+1,  -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v2,    -area * dy));
        triplets.push_back(Eigen::Triplet<double>(2*v2+1, 2*v2+1,  -area * (-dx)));
    }
}

/**
 * Example usage
 */
void example_usage() {
    // Suppose we have an island with 4 vertices
    int num_vertices = 4;

    // Triangle with local indices 0, 1, 2
    int v0 = 0, v1 = 1, v2 = 2;

    // 3D positions (example)
    Vec3 p0(0.0, 0.0, 0.0);
    Vec3 p1(1.0, 0.0, 0.0);
    Vec3 p2(0.5, 0.866, 0.0);  // Equilateral triangle

    // Build matrix
    std::vector<Eigen::Triplet<double>> triplets;

    add_triangle_lscm_contribution(triplets, v0, v1, v2, p0, p1, p2);

    // Create sparse matrix (2n × 2n)
    Eigen::SparseMatrix<double> A(2 * num_vertices, 2 * num_vertices);
    A.setFromTriplets(triplets.begin(), triplets.end());

    std::cout << "Matrix size: " << A.rows() << " × " << A.cols() << std::endl;
    std::cout << "Non-zero entries: " << A.nonZeros() << std::endl;

    // To complete LSCM:
    // 1. Loop over all triangles, calling add_triangle_lscm_contribution
    // 2. Pin 2 boundary vertices
    // 3. Solve the system
}

int main() {
    example_usage();
    return 0;
}
