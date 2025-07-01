#include "recursive.h"
#include "draw.h"
#include "bvh_interface.h"
#include "intersect.h"
#include "extra.h"
#include "light.h"
#include <iostream>
#include "bvh.h"
#include <glm/gtx/string_cast.hpp>



// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
// This function is provided as-is. You do not have to implement it.
// Given a range of rays, render out all rays and average the result
glm::vec3 renderRays(RenderState& state, std::span<const Ray> rays, int rayDepth)
{
    // 1) If DOF is not enabled, just accumulate color
    if (!state.features.extra.enableDepthOfField) {
        glm::vec3 L(0.f);
        for (const auto& ray : rays) {
            L += renderRay(state, ray, rayDepth);
        }
        return L / float(rays.size());
    }

    // 2) DOF is enabled: accumulate color
    glm::vec3 L(0.f);

    // 3) If debug drawing is off, skip geometry draws
    if (!state.features.enableDebugDraw) {
        for (const auto& ray : rays) {
            L += renderRay(state, ray, rayDepth);
        }
        return L / float(rays.size());
    }

    // ---------------------------------------------------------------------
    // 4) Retrieve DOF parameters & compute average lens/focal/image centers
    float lensRadius         = state.features.extra.dofSettings.lendRadius;
    float focalDistance      = state.features.extra.dofSettings.focalDistance;
    float imagePlaneDistance = state.features.extra.dofSettings.imagePlaneDistance;

    glm::vec3 lensPlaneCenter(0.f);
    if (!debugData.lensPoints.empty()) {
        for (auto& p : debugData.lensPoints) {
            lensPlaneCenter += p;
        }
        lensPlaneCenter /= float(debugData.lensPoints.size());
    }

    glm::vec3 focalPlaneCenter(0.f);
    if (!debugData.focalPoints.empty()) {
        for (auto& p : debugData.focalPoints) {
            focalPlaneCenter += p;
        }
        focalPlaneCenter /= float(debugData.focalPoints.size());
    }

    glm::vec3 imagePlaneCenter(0.f);
    if (!debugData.imagePlanePoints.empty()) {
        for (auto& p : debugData.imagePlanePoints) {
            imagePlaneCenter += p;
        }
        imagePlaneCenter /= float(debugData.imagePlanePoints.size());
    } else {
        // Fallback if no image-plane points
        imagePlaneCenter = lensPlaneCenter + imagePlaneDistance * debugData.cameraForward;
    }

    // ---------------------------------------------------------------------
    // (A) Check that all focal points lie in one plane 
    if (debugData.focalPoints.size() >= 3) {
        glm::vec3 p0 = debugData.focalPoints[0];
        glm::vec3 p1 = debugData.focalPoints[1];
        glm::vec3 p2 = debugData.focalPoints[2];
        glm::vec3 planeNormal = glm::normalize(glm::cross(p1 - p0, p2 - p0));

        float maxDeviation = 0.0f;
        for (size_t i = 3; i < debugData.focalPoints.size(); i++) {
            glm::vec3 v = debugData.focalPoints[i] - p0;
            float dist  = glm::dot(v, planeNormal);
            maxDeviation = glm::max(maxDeviation, glm::abs(dist));
        }

    }

    // ---------------------------------------------------------------------
    // (C) Draw debug geometry
    //    - a thin lens (red) for the lens plane,
    //    - a green plane for the focal plane,
    //    - a blue plane for the image plane.

    // 1) Draw the lens as a circle (using drawThinLens)
    drawThinLens(
        lensPlaneCenter,
        debugData.cameraForward, // the lens normal
        lensRadius,  // radiusX
        lensRadius,  // radiusY
        glm::vec3(1, 0, 0) // red
    );

    // 2) Draw the focal plane in green
    //    We'll just call drawPlane(...) with the plane center, normal=forward, up=debugData.cameraUp
    drawPlane(
        focalPlaneCenter,
        debugData.cameraForward,
        debugData.cameraUp,
        glm::vec3(0, 1, 0) // green
    );

    // 3) Draw the image plane in blue
    drawPlane(
        imagePlaneCenter,
        debugData.cameraForward,
        debugData.cameraUp,
        glm::vec3(0, 0, 1) // blue
    );

    // ---------------------------------------------------------------------
    // (D) Trace each ray for shading, also draw debug lines/spheres
    for (size_t i = 0; i < rays.size(); i++) {
        // accumulate color
        L += renderRay(state, rays[i], rayDepth);

        // lens point: small red sphere
        if (i < debugData.lensPoints.size()) {
            drawSphere(debugData.lensPoints[i], 0.01f, glm::vec3(1, 0, 0));
        }
        // focal point: yellow sphere
        if (i < debugData.focalPoints.size()) {
            drawSphere(debugData.focalPoints[i], 0.02f, glm::vec3(1, 1, 0));
        }

        // Ray from image-plane point to lens (blue line)
        if (i < debugData.imagePlanePoints.size() && i < debugData.lensPoints.size()) {
            glm::vec3 start = debugData.imagePlanePoints[i];
            glm::vec3 end   = debugData.lensPoints[i];
            glm::vec3 dir   = glm::normalize(end - start);
            float dist      = glm::length(end - start);
            drawRay(Ray(start, dir * dist), glm::vec3(0, 0, 1));
        }

        // Ray from lens to focal point (yellow line)
        if (i < debugData.lensPoints.size() && i < debugData.focalPoints.size()) {
            glm::vec3 start = debugData.lensPoints[i];
            glm::vec3 end   = debugData.focalPoints[i];
            glm::vec3 dir   = glm::normalize(end - start);
            float dist      = glm::length(end - start);
            drawRay(Ray(start, dir * dist), glm::vec3(1, 1, 0));
        }
    }
 
    // ---------------------------------------------------------------------
    // 6) Return the averaged color from all DOF rays
    return L / float(rays.size());
}


// This method is provided as-is. You do not have to implement it.
// Given a camera ray (or secondary ray), tests for a scene intersection, and
// dependent on the results, evaluates the following functions which you must
// implement yourself:
// - `computeLightContribution()` and its submethods
// - `renderRaySpecularComponent()`, `renderRayTransparentComponent()`, `renderRayGlossyComponent()`
glm::vec3 renderRay(RenderState& state, Ray ray, int rayDepth)
{
    // Trace the ray into the scene. If nothing was hit, return early
    HitInfo hitInfo;
    if (!state.bvh.intersect(state, ray, hitInfo)) {
        if (state.features.enableDebugDraw) {
            drawRay(ray, glm::vec3(1, 0, 0));
        }
        return sampleEnvironmentMap(state, ray);
    }

    // Return value: the light along the ray
    // Given an intersection, estimate the contribution of scene lights at this intersection
    glm::vec3 Lo = computeLightContribution(state, ray, hitInfo);

    // Draw an example debug ray for the incident ray (feel free to modify this for yourself)
    drawRay(ray, glm::vec3(1.0f));

    // Given that recursive components are enabled, and we have not exceeded maximum depth,
    // estimate the contribution along these components
    if (rayDepth < 6) {
        bool isReflective = glm::any(glm::notEqual(hitInfo.material.ks, glm::vec3(0.0f)));
        bool isTransparent = hitInfo.material.transparency != 1.f;

        // Default, specular reflections
        if (state.features.enableReflections && !state.features.extra.enableGlossyReflection && isReflective) {
            renderRaySpecularComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Alternative, glossy reflections
        if (state.features.enableReflections && state.features.extra.enableGlossyReflection && isReflective) {
            renderRayGlossyComponent(state, ray, hitInfo, Lo, rayDepth);
        }

        // Transparency passthrough
        if (state.features.enableTransparency && isTransparent) {
            renderRayTransparentComponent(state, ray, hitInfo, Lo, rayDepth);
        }
    }

    return Lo;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a mirrored ray
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a reflected ray
// This method is unit-tested, so do not change the function signature.
Ray generateReflectionRay(Ray ray, HitInfo hitInfo)
{
    // TODO: generate a mirrored ray
    //       if you use glm::reflect, you will not get points for this method!
    Ray r {};

    r.origin = ray.origin+ray.direction*ray.t;
    r.direction = -hitInfo.normal * glm::dot(hitInfo.normal, ray.direction);
    r.direction *= 2.0;
    r.direction += ray.direction;
    
    r.t = FLT_MAX;
    r.direction = glm::normalize(r.direction);
    r.origin += r.direction * 0.00001f;
    return r;
}

// TODO: Standard feature
// Given an incident ray and a intersection point, generate a passthrough ray for transparency,
// starting at the intersection point and continuing in the same direction.
// - Ray;     the indicent ray
// - HitInfo; hit struct for the intersection point
// - return;  a passthrough ray for transparency
// This method is unit-tested, so do not change the function signature.
Ray generatePassthroughRay(Ray ray, HitInfo hitInfo)
{

    // TODO: generate a passthrough ray
    Ray r {};
    r.direction = ray.direction;
    r.t = FLT_MAX;
    r.origin = ray.origin + ray.direction * (0.00001f+ray.t);
    return r;
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a mirrored ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and adding the result times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRaySpecularComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generateReflectionRay()
    Ray r = generateReflectionRay(ray, hitInfo);
    
    glm::vec3 col = renderRay(state, r, rayDepth+1);
    col = glm::clamp(col, { 0, 0, 0 }, { 1, 1, 1 });
    hitColor += col * hitInfo.material.ks;
    if (state.features.enableDebugDraw) {
        glm::vec3 secondaryColour = col * float(std::pow(0.9f, rayDepth+1)) + glm::vec3 {0.1f, 0.1f, 0.1f};
        HitInfo newHit;
        state.bvh.intersect(state, r, newHit);
        drawRay(r, secondaryColour);

        Ray c {};
        c.origin = r.origin;
        c.direction = hitInfo.normal;
        c.t = 1;

        //normal
        drawRay(c, glm::vec3{1, 0, 1});
    }
    
    // ...
}

// TODO: standard feature
// Given a camera ray (or secondary ray) and an intersection, evaluates the contribution
// of a passthrough transparent ray, recursively evaluating renderRay(..., depth + 1) along this ray,
// and correctly alpha blending the result with the current intersection's hit color
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is unit-tested, so do not change the function signature.
void renderRayTransparentComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // TODO; you should first implement generatePassthroughRay()
    Ray r = generatePassthroughRay(ray, hitInfo);
    glm::vec3 col = renderRay(state, r, rayDepth + 1);
    HitInfo info;

    state.bvh.intersect(state, r, info);

    hitColor = (1 - hitInfo.material.transparency) * col + (hitInfo.material.transparency) * hitColor; 

    if (state.features.enableDebugDraw && !state.features.enableShadows) {
        drawRay(r, hitColor);

        /*
        I chose to illustrate the transparency by the length of the normal coming from the face. 
        The length varies from 1 to 0 depending on the transparency, the more transparent the face 
        is, the longer the normal. I am using the same normal to demonstrate the kd of the material.
        */

        Ray vis;
        vis.origin = ray.origin + ray.direction * ray.t;
        vis.direction = hitInfo.normal;
        vis.t = hitInfo.material.transparency;

        drawRay(vis, hitInfo.material.kd);
    }
 
    // ...
}