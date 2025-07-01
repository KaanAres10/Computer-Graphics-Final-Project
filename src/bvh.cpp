#pragma once
#include "bvh.h"
#include <stack>
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "interpolate.h"
#include <vector>
#include "draw.h"

#include "texture.h"

glm::vec3 hsvToRgb(float h, float s, float v)
{
    float c = v * s; // Chroma
    float x = c * (1.0f - fabs(fmod(h / 60.0f, 2) - 1.0f));
    float m = v - c;

    float r, g, b;
    if (h >= 0 && h < 60) {
        r = c;
        g = x;
        b = 0;
    } else if (h >= 60 && h < 120) {
        r = x;
        g = c;
        b = 0;
    } else if (h >= 120 && h < 180) {
        r = 0;
        g = c;
        b = x;
    } else if (h >= 180 && h < 240) {
        r = 0;
        g = x;
        b = c;
    } else if (h >= 240 && h < 300) {
        r = x;
        g = 0;
        b = c;
    } else {
        r = c;
        g = 0;
        b = x;
    }

    return glm::vec3(r + m, g + m, b + m);
}

// TODO Standard Feature
// Hierarchy traversal routine; you must implement this method and implement it carefully!
//
// The default implementation uses precompiled methods and only performs rudimentary updates to hitInfo.
// For correct normal interpolation, barycentric coordinates, etc., you have to implement this traversal yourself.
// If you have implemented the traversal method for assignment 4B, you might be able to reuse most of it.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
// - state;    current render state (containing scene, features, ...)
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
bool BVH::intersect(RenderState& state, Ray& ray, HitInfo& hitInfo) const
{
    if (!state.features.enableAccelStructure) {
        return intersectRayWithBVH(state, *this, ray, hitInfo);
    }
    if (state.features.enableDebugDraw) {
        ray.t = FLT_MAX;
        drawRay(ray);
    }

    glm::vec2 hitPos;


    // TODO: implement here your (probably stack-based) BVH traversal.
    // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
    // data is not easily extracted. Helper methods are available, however:
    // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
    // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
    // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
    //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
    //
    // In short, you will have to step down the bvh, node by node, and intersect your ray
    // with the node's AABB. If this intersection passes, you should:
    // - if the node is a leaf, intersect with the leaf's primitives
    // - if the node is not a leaf, test the left and right children as well!
    //
    // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
    // and it is likewise possible for a ray to hit both children of a node.
    //
    // Make sure to update the hitInfo in case of a hit
    std::stack<std::pair<int, BVHInterface::Node>> stack;
    if (!intersectRayWithShape(m_nodes[0].aabb, ray)) {
        ray.t = FLT_MAX;
        return false;
    }
        
    stack.push({ 0, m_nodes[0] });
    bool hitTriangle = false;
    float rayDistance = FLT_MAX;
    std::vector<std::pair<int, Node>> visited;
    while (!stack.empty()) {
        BVHInterface::Node node = stack.top().second;
        int level = stack.top().first;
        stack.pop();
        if (state.features.enableDebugDraw)
            visited.push_back({ level, node });

        if (node.isLeaf()) {
            for (int i = 0; i < node.primitiveCount(); i++) {
                Primitive triangle = m_primitives[node.primitiveOffset() + i];
                
                ray.t = FLT_MAX;
                if (!intersectRayWithTriangle(triangle.v0.position, triangle.v1.position, triangle.v2.position, ray, hitInfo)) {
                    continue;
                }
                    
                if (state.features.enableDebugDraw) {
                    drawTriangle(triangle.v0, triangle.v1, triangle.v2, { 1.0, 0.1, 1.0 });
                    drawSphere(ray.origin + ray.direction * ray.t, 0.005f, {1.0f, 1.0f, 1.0f});
                }
                if (rayDistance > ray.t) {
                    hitTriangle = true;
                    hitInfo.material = state.scene.meshes[triangle.meshID].material;
                    rayDistance = ray.t;
                    hitInfo.barycentricCoord = computeBarycentricCoord(triangle.v0.position, triangle.v1.position, triangle.v2.position, ray.origin + ray.direction * ray.t);
                    hitInfo.texCoord = interpolateTexCoord(triangle.v0.texCoord, triangle.v1.texCoord, triangle.v2.texCoord, hitInfo.barycentricCoord);
                    // Interpolate and normalize the normal
                    glm::vec3 interpolatedNormal = interpolateNormal(
                        triangle.v0.normal,
                        triangle.v1.normal,
                        triangle.v2.normal,
                        hitInfo.barycentricCoord);

                    // Flip the normal if it's back-facing
                    if (glm::dot(ray.direction, interpolatedNormal) > 0.0f) {
                        if (state.features.enableShadows) {
                            interpolatedNormal = -interpolatedNormal;
                        }
                    }
                    hitInfo.normal = glm::normalize(interpolatedNormal);

                    hitPos = glm::vec2(hitInfo.texCoord.x, 1.0f - hitInfo.texCoord.y);
                }
            }
            continue;
        }

        ray.t = FLT_MAX;
        if (intersectRayWithShape(m_nodes[node.leftChild()].aabb, ray))
            stack.push({ level + 1, m_nodes[node.leftChild()] });

        ray.t = FLT_MAX;
        if (intersectRayWithShape(m_nodes[node.rightChild()].aabb, ray))
            stack.push({ level + 1, m_nodes[node.rightChild()] });
    }

    for (int i = 0; i < visited.size(); i++) {
        float h = float(i) / ((float)visited.size() - 1) * 180;
        float depth = float(visited[i].first) / (float)m_numLevels * 0.9 + 0.1;
        drawAABB(visited[i].second.aabb, DrawMode::Wireframe, hsvToRgb(h, 1, depth));
    }

    if (state.features.enableTextureMapping && state.features.enableDebugDraw && hitTriangle && hitInfo.material.kdTexture != nullptr) {
        textureVisualDebug(hitInfo.material.kdTexture, hitPos, state);
    }

    ray.t = rayDistance;
    return hitTriangle;
}
