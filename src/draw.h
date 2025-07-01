#pragma once
#include "scene.h"
#include <framework/mesh.h>
#include <framework/ray.h>
#include <utility> // std::forward

// Add your own custom visual debug draw functions here then implement it in draw.cpp.
// You are free to modify the example one however you like.
void drawExampleOfCustomVisualDebug();

void drawRay(const Ray& ray, const glm::vec3& color = glm::vec3(1.0f));

void drawAABB(const AxisAlignedBox& box, DrawMode drawMode = DrawMode::Filled, const glm::vec3& color = glm::vec3(1.0f), float transparency = 1.0f);

void drawTriangle (const Vertex& v0, const Vertex& v1, const Vertex& v2 );
void drawMesh(const Mesh& mesh);
void drawSphere(const Sphere& sphere);
void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color = glm::vec3(1.0f));
void drawScene(const Scene& scene);

void drawPlane(const glm::vec3& position, const glm::vec3& normal, const glm::vec3& up, const glm::vec3& color);

void drawEllipse(const glm::vec3& position, const glm::vec3& normal, const glm::vec3& up, float radiusX, float radiusY, const glm::vec3& color);

void drawLens(const glm::vec3& position, const glm::vec3& normal, float radius, const glm::vec3& color);

void drawThinLens(const glm::vec3& position, const glm::vec3& normal, float radiusX, float radiusY, const glm::vec3& color);
void drawTriangle(const Vertex& v0, const Vertex& v1, const Vertex& v2, const glm::vec3& colour);
void drawLine(const glm::vec3& pos, const glm::vec3& pos2);
void drayRayYawPitch(const glm::vec3& origin, const glm::vec2& rotation, float length = 0.1f);
void drawSpline(const BSpline& spline, float step=0.01f);
