#include "render.h"
#include "texture.h"
#include <cmath>
#include <fmt/core.h>
#include <glm/geometric.hpp>
#include <glm/gtx/string_cast.hpp>
#include <shading.h>
#include <ranges>
#include <algorithm>
#include <iostream>

// This function is provided as-is. You do not have to implement it (unless
// you need to for some extra feature).
// Given render state and an intersection, based on render settings, sample
// the underlying material data in the expected manner.
glm::vec3 sampleMaterialKd(RenderState& state, const HitInfo& hitInfo)
{
    if (state.features.enableTextureMapping && hitInfo.material.kdTexture) {
        if (state.features.enableBilinearTextureFiltering) {
            return sampleTextureBilinear(*hitInfo.material.kdTexture, hitInfo.texCoord);
        } else {
            return sampleTextureNearest(*hitInfo.material.kdTexture, hitInfo.texCoord);
        }
    } else {
        return hitInfo.material.kd;
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate the scene-selected shading model, returning the reflected light towards the target.
glm::vec3 computeShading(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo)
{
    // Hardcoded color gradient. Feel free to modify this
    static LinearGradient gradient = {
        .components = {
            {-1.0f, glm::vec3(1.0f, 0.0f, 0.0f) }, // Red - Phong Stronger
            { 1.0f, glm::vec3(0.0f, 1.0f, 1.0f) } // Cyan - Blinn-Phong Stronger
        }
    };

    if (state.features.enableShading) {
        const glm::vec3 kd = sampleMaterialKd(state, hitInfo);
        switch (state.features.shadingModel) {
            case ShadingModel::Lambertian:
                return computeLambertianModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
            case ShadingModel::Phong:
                return computePhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
            case ShadingModel::BlinnPhong:
                return computeBlinnPhongModel(state, cameraDirection, lightDirection, lightColor, hitInfo, kd);
            case ShadingModel::LinearGradient:
                return computeLinearGradientModel(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
            case ShadingModel::LinearGradientComparison:
                return computeLinearGradientModelComparison(state, cameraDirection, lightDirection, lightColor, hitInfo, gradient);
            };
    }

    return lightColor * sampleMaterialKd(state, hitInfo);
}

// TODO: Standard feature
// Given a number ti between [-1, 1], sample from the gradient's components and return the
// linearly interpolated color, for which ti lies in the interval between the t-values of two
// components, or on a boundary. If ti falls outside the gradient's smallest/largest components,
// the nearest component must be sampled.
// - ti; a number between [-1, 1]
// This method is unit-tested, so do not change the function signature.
glm::vec3 LinearGradient::sample(float ti) const
{
    if (components.empty())
        return glm::vec3(0.5f); // If no components exist return gray

    // Sort components by t value
    std::vector<Component> sortedComponents = components;
    std::sort(sortedComponents.begin(), sortedComponents.end(), [](const Component& a, const Component& b) {
        return a.t < b.t;
    });

    // Clamp to [-1,1]
    ti = glm::clamp(ti, -1.0f, 1.0f);

    // If ti is smaller than smallest component return smallest component of boundary
    if (ti <= sortedComponents.front().t)
        return sortedComponents.front().color;

    // If ti is larger than largest component return largest component of boundary
    if (ti >= sortedComponents.back().t)
        return sortedComponents.back().color;

    for (size_t i = 0; i < sortedComponents.size() - 1; i++)
    {
        Component& c1 = sortedComponents[i];
        Component& c2 = sortedComponents[i + 1];

        if (ti >= c1.t && ti <= c2.t)
        {
            // Linear interpolation 
            float alpha = (ti - c1.t) / (c2.t - c1.t);
            return (1.0f - alpha) * c1.color + alpha * c2.color;
        }
    }

    return glm::vec3(0.0f); // Should never reach
}

// TODO: Standard feature
// Given a camera direction, a light direction, a relevant intersection, and a color coming in
// from the light, evaluate a diffuse shading model, such that the diffuse component is sampled not
// from the intersected material, but a provided linear gradient, based on the cosine of theta
// as defined in the diffuse shading part of the Phong model.
//
// - state;           the active scene, feature config, and the bvh
// - cameraDirection; exitant vector towards the camera (or secondary position)
// - lightDirection;  exitant vector towards the light
// - lightColor;      the color of light along the lightDirection vector
// - hitInfo;         hit object describing the intersection point
// - gradient;        the linear gradient object
// - return;          the result of shading
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeLinearGradientModel(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    float cos_theta = glm::dot(lightDirection, hitInfo.normal);
    cos_theta = glm::clamp(cos_theta, -1.0f, 1.0f);
    glm::vec3 gradient_color = gradient.sample(cos_theta);
    return gradient_color * lightColor;
}

glm::vec3 computeLinearGradientModelComparison(RenderState& state, const glm::vec3& cameraDirection, const glm::vec3& lightDirection, const glm::vec3& lightColor, const HitInfo& hitInfo, const LinearGradient& gradient)
{
    // Phong
    glm::vec3 reflectDir = glm::normalize(-lightDirection + 2 * glm::dot(lightDirection, hitInfo.normal) * (hitInfo.normal));
    glm::vec3 phongSpecular = lightColor * hitInfo.material.ks * 
                              glm::pow(glm::max(glm::dot(reflectDir, glm::normalize(cameraDirection)), 0.0f), hitInfo.material.shininess);


    // Blinn Phong
    glm::vec3 halfWayDir = glm::normalize(lightDirection + cameraDirection);
    glm::vec3 blinnPhongSpecular = lightColor * hitInfo.material.ks * 
                                    glm::pow(glm::max(glm::dot(halfWayDir, hitInfo.normal), 0.0f), hitInfo.material.shininess); 

    // Differnce
    float difference = glm::length(phongSpecular - blinnPhongSpecular);
    float maxSpecular = glm::max(glm::length(phongSpecular), glm::length(blinnPhongSpecular) + 1e-5f);
    difference = glm::clamp(difference / maxSpecular, 0.0f, 1.0f);

    float ti;
    if (glm::length(phongSpecular) > glm::length(blinnPhongSpecular))
    {
        ti = -difference;
    }
    else
    {
        ti = difference;
    }
    ti = glm::sign(ti) * glm::pow(glm::abs(ti), 1.5f);

    return gradient.sample(ti);
}
