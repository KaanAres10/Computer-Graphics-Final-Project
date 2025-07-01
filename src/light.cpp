#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <iostream>
DISABLE_WARNINGS_POP()


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{
    position = light.endpoint0 + sample * (light.endpoint1 - light.endpoint0);
    color = light.color0 + sample * (light.color1 - light.color0);
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    // TODO: implement this function.
    position = light.v0 + sample.x * light.edge01 + sample.y * light.edge02;
    
    glm::vec3 c0 = light.color0;
    glm::vec3 c1 = light.color1;
    glm::vec3 c2 = light.color2;
    glm::vec3 c3 = light.color3;

    glm::vec3 colorEdge1_x = glm::mix(c0, c1, sample.x);
    glm::vec3 colorEdge2_x = glm::mix(c2, c3, sample.x);

    color = glm::mix(colorEdge1_x, colorEdge2_x, sample.y);
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer
        // TODO: implement this function; currently, the light simply passes through
        glm::vec3 intersectionPoint = ray.origin + ray.t * (ray.direction);
        glm::vec3 lightDir = glm::normalize(lightPosition - intersectionPoint);

        // Prevent self shadowing
        glm::vec3 offSetOrigin = intersectionPoint + hitInfo.normal * 1e-3f;

        Ray shadowRay(offSetOrigin, lightDir);

        HitInfo shadowHit;

        float epsilon = 1e-4f;



        if (state.bvh.intersect(state,shadowRay, shadowHit))
        {
            glm::vec3 shadowHitLocation = shadowRay.origin + shadowRay.t * (shadowRay.direction);
            float distToBlock = glm::length(shadowHitLocation - offSetOrigin);
            float distToLight = glm::length(lightPosition - offSetOrigin);  
            
            //// Debug
            if (state.features.enableDebugDraw)
            {

                if (distToBlock < distToLight - epsilon)
                {
                    drawRay(shadowRay, lightColor);
                    drawSphere(shadowHitLocation, 0.01f, glm::vec3(1.0f, 0.0f, 0.0f)); // Draw blocking point
                } 

            }

            if (distToBlock < distToLight - epsilon)
            {
                return false;
            }
        }
        
        return true;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: implement this function; currently, the light simply passes through

    glm::vec3 lightPassed = lightColor;

    if (!state.features.enableShadows)
        return lightPassed;

    glm::vec3 intersectionPoint = ray.origin + ray.t * (ray.direction);
    glm::vec3 lightDir = glm::normalize(lightPosition - intersectionPoint);


    // Prevent self shadowing
    glm::vec3 offSetOrigin = intersectionPoint + hitInfo.normal * 1e-3f;

    Ray shadowRay(offSetOrigin, lightDir);
    HitInfo shadowHit;

    if (state.bvh.intersect(state, shadowRay, shadowHit)) {
        glm::vec3 shadowHitLocation = shadowRay.origin + shadowRay.t * (shadowRay.direction);
        float distToBlock = glm::length(shadowHitLocation - offSetOrigin);
        float distToLight = glm::length(lightPosition - offSetOrigin);


        if (distToBlock < distToLight) {
            float alpha = shadowHit.material.transparency;
            lightPassed = lightPassed * shadowHit.material.kd * (1.0f - alpha);

            // Correct Version
            //lightPassed = lightPassed * shadowHit.material.kd * (1.0f - alpha) + alpha * lightPassed;
           
            // Debug Visually Ray attenuation 
            if (state.features.enableDebugDraw) {
                if (glm::length(lightPassed) > 0 && alpha > 0.0f)
                {
                    drawRay(Ray(shadowRay.origin, glm::normalize(shadowHitLocation - shadowRay.origin)), lightPassed);
                    drawSphere(shadowHitLocation, 0.01f, lightPassed);
                }
            }

            if (alpha == 0.0f)
            {
                lightPassed = glm::vec3(0.0f);
            }
            shadowRay.origin = shadowHitLocation + shadowHit.normal * 1e-3f;
        }
        else
        {
        }
    }

    // Debug Visually Final Color
    if (state.features.enableDebugDraw) {
        if (glm::length(lightPassed) > 0)
            drawRay(Ray(shadowRay.origin, lightDir), lightPassed);
    }

    return lightPassed;
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    // TODO: modify this function to incorporate visibility corerctly
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;

    glm::vec3 visibility = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);

    // Debugging
    if (state.features.enableDebugDraw)
    {
        glm::vec3 rayColor = glm::length(visibility) > 0 ? glm::vec3(0.0f, 1.0f, 0.0f) : glm::vec3(1.0f, 0.0f, 0.0f); // If visible return green else return red
        drawRay(Ray(p, glm::normalize(light.position - p)), rayColor); 
    }

    if (glm::length(visibility) > 0)
    {
       return visibility * computeShading(state, v, l, light.color, hitInfo); 
    } 

    return glm::vec3(0);
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model

    glm::vec3 accumulatedLight = glm::vec3(0);
    
    glm::vec3 intersectionPoint = ray.origin + ray.t * (ray.direction);


    for (uint32_t i = 0; i < numSamples; i++)
    {

        // Sample
        glm::vec3 lightPos, lightColor;
        sampleSegmentLight(state.sampler.next_1d(), light, lightPos, lightColor);


        // Debug Visually
        if (state.features.enableDebugDraw)
        {
            drawSphere(lightPos, 0.005f, lightColor);
        }
       
        // Visibility
        glm::vec3 visibility = visibilityOfLightSample(state, lightPos, lightColor, ray, hitInfo);


        if (state.features.enableDebugDraw)
        {
            if (glm::length(visibility) > 0)
            {
                drawRay(Ray(intersectionPoint, glm::normalize(lightPos - intersectionPoint)), visibility);
            }

        }


        // If visible
        if (glm::length(visibility) > 0)
        {
            glm::vec3 lightDir = glm::normalize(lightPos - intersectionPoint);
            glm::vec3 viewDir = -ray.direction;
            accumulatedLight += visibility * computeShading(state, viewDir, lightDir, visibility, hitInfo);
        }

    }


    return accumulatedLight / static_cast<float>(numSamples);
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model

    bool testDiff = false;

    /*
      Compare samples drawn with the provided sampler for 2D to the samples that you obtain by choosing both
      coordinates of a sample with the 1D random sampling
    */
    if (testDiff) {
        for (uint32_t i = 0; i < numSamples; i++) {

            // Sample 2D together
            glm::vec3 lightPos2D, lightColor2D;
            sampleParallelogramLight(state.sampler.next_2d(), light, lightPos2D, lightColor2D);

            drawSphere(lightPos2D, 0.005f, glm::vec3(0.0f, 1.0f, 0.0f)); // Green 2D together


            // Sample 1D seperate 
            glm::vec3 lightPos1D, lightColor1D;
            float x = state.sampler.next_1d();
            float y = state.sampler.next_1d();

            // 1d for x, 1d for y, seperate
            sampleParallelogramLight({ x, y }, light, lightPos1D, lightColor1D);
            drawSphere(lightPos1D, 0.005f, glm::vec3(1.0f, 0.0f, 0.0f)); // Red 1D together

        }
        return glm::vec3(0);
    } 

    glm::vec3 accumulatedLight = glm::vec3(0);

    glm::vec3 intersectionPoint = ray.origin + ray.t * (ray.direction);

    for (uint32_t i = 0; i < numSamples; i++) {

        // Sample
        glm::vec3 lightPos, lightColor;

        float x = state.sampler.next_1d();
        float y = state.sampler.next_1d();


        // 1d for x, 1d for y, seperate
        //sampleParallelogramLight({ x, y }, light, lightPos, lightColor);

        // 2d together
        sampleParallelogramLight(state.sampler.next_2d(), light, lightPos, lightColor);

        // Debug Visually
        if (state.features.enableDebugDraw) {
            drawSphere(lightPos, 0.005f, lightColor);
        }

        // Visibility
        glm::vec3 visibility = visibilityOfLightSample(state, lightPos, lightColor, ray, hitInfo);

        if (state.features.enableDebugDraw) {
            if (glm::length(visibility) > 0) {
                drawRay(Ray(intersectionPoint, glm::normalize(lightPos - intersectionPoint)), visibility);
            }
        }

        // If visible
        if (glm::length(visibility) > 0) {
            glm::vec3 lightDir = glm::normalize(lightPos - intersectionPoint);
            glm::vec3 viewDir = -ray.direction;
            accumulatedLight +=  visibility * computeShading(state, viewDir, lightDir, visibility, hitInfo);
        }
    }

    return accumulatedLight / static_cast<float>(numSamples);
}


// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }

    return Lo;
}