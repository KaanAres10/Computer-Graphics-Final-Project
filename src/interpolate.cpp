#include "interpolate.h"
#include <glm/geometric.hpp>

// TODO Standard feature
// Given three triangle vertices and a point on the triangle, compute the corresponding barycentric coordinates of the point.
// and return a vec3 with the barycentric coordinates (alpha, beta, gamma).
// - v0;     Triangle vertex 0
// - v1;     Triangle vertex 1
// - v2;     Triangle vertex 2
// - p;      Point on triangle
// - return; Corresponding barycentric coordinates for point p.
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeBarycentricCoord(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& p)
{
    // TODO: implement this function.

    glm::vec3 side1 = v1 - v0;
    glm::vec3 side2 = v2 - v0;

    float area = glm::length(glm::cross(side1, side2));
    float areaA = glm::length(glm::cross(v1 - p, v2 - p));
    float areaB = glm::length(glm::cross(v2 - p, v0 - p));
    float areaC = glm::length(glm::cross(v0 - p, v1 - p));

    float a = areaA / area;
    float b = areaB / area;
    float c = areaC / area;

    return glm::vec3 { a, b, c };
}

// TODO Standard feature
// Linearly interpolate three normals using barycentric coordinates.
// - n0;     Triangle normal 0
// - n1;     Triangle normal 1
// - n2;     Triangle normal 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated normal.
// This method is unit-tested, so do not change the function signature.
glm::vec3 interpolateNormal(const glm::vec3& n0, const glm::vec3& n1, const glm::vec3& n2, const glm::vec3 bc)
{
    // TODO: implement this function.
    glm::vec3 interp = n0 * bc.x + n1 * bc.y + n2 * bc.z;
    return glm::normalize(interp);
}

// TODO Standard feature
// Linearly interpolate three texture coordinates using barycentric coordinates.
// - n0;     Triangle texture coordinate 0
// - n1;     Triangle texture coordinate 1
// - n2;     Triangle texture coordinate 2
// - bc;     Barycentric coordinate
// - return; The smoothly interpolated texturre coordinate.
// This method is unit-tested, so do not change the function signature.
glm::vec2 interpolateTexCoord(const glm::vec2& t0, const glm::vec2& t1, const glm::vec2& t2, const glm::vec3 bc)
{
// TODO: implement this function.

    return bc.x * t0 + bc.y * t1 + bc.z * t2;
}
