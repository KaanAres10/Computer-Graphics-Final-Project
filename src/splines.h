#pragma once
#include <vector>
#include <glm/ext/vector_float3.hpp>
#include <glm/ext/vector_float2.hpp>

enum SplineType {
    Linear = 0,
    Quadratic = 1,
    Cubic = 2,
};

struct SplinePoint {
    glm::vec3 position;
    glm::vec2 rotation;
};

struct BSpline {
    std::vector<glm::vec3> position;
    std::vector<glm::vec2> rotation;
    bool repeatKnots;
    SplineType type;
    float t;
    int rayNumber;
};

SplinePoint interpolateBSpline(const BSpline& bSpline, float t);

template <typename T>
T interpolateUniformLinear(const std::vector<T>& controlPoints, float t);

template <typename T>
T interpolateUniformQuadratic(const std::vector<T>& controlPoints, float t);

template <typename T>
T interpolateUniformCubic(const std::vector<T>& controlPoints, float t);

template <typename T>
T interpolateConnectedQuadratic(const std::vector<T>& controlPoints, float t);

template <typename T>
T interpolateConnectedCubic(const std::vector<T>& controlPoints, float t);