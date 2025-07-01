#include "extra.h"
#include "bvh.h"
#include "light.h"
#include "recursive.h"
#include "shading.h"
#include "render.h"
#include <framework/trackball.h>
#include <random>
#include <draw.h>
#include <screen.cpp>
#include <glm/gtx/string_cast.hpp>

DebugRayData debugData; // Define the global variable


glm::vec2 sampleDisk(float radius) {
    float r = radius * std::sqrt(static_cast<float>(rand()) / RAND_MAX);
    float theta = 2.0f * glm::pi<float>() * (static_cast<float>(rand()) / RAND_MAX);
    return glm::vec2(r * std::cos(theta), r * std::sin(theta));
}


// Generare depth of field ray that focuses on focal point
// 

#include <imgui/imgui.h>
#include <GLFW/glfw3.h>
#include "draw.h"
#include "texture.h"
#include "intersect.h"
#include "interpolate.h"

// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of Depth of Field. Here, you generate camera rays s.t. a focus point and a thin lens camera model
// are in play, allowing objects to be in and out of focus.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithDepthOfField(
    const Scene& scene,
    const BVHInterface& bvh,
    const Features& features,
    const Trackball& camera,
    Screen& screen)
{
    // Only do DOF if the flag is enabled
    if (!features.extra.enableDepthOfField) {
        return;
    }

    auto& pixels = screen.pixels();
    glm::ivec2 resolution = screen.resolution();

    // Number of lens samples per pixel
    const int samples = features.numPixelSamples;

    // For each pixel
    for (int y = 0; y < resolution.y; y++) {
        for (int x = 0; x < resolution.x; x++) {

            // Convert pixel coords to NDC in range [-1..1]
            float ndcX = ((x + 0.5f) / float(resolution.x)) * 2.0f - 1.0f;
            float ndcY = 1.0f - ((y + 0.5f) / float(resolution.y)) * 2.0f;

            glm::vec2 pixelCenterNDC(ndcX, ndcY);

            // Generate multiple DOF rays through this pixel
            std::vector<Ray> dofRays = generateDepthOfFieldRays(
                camera,
                pixelCenterNDC,
                resolution,
                features,
                samples);

            // Trace/accumulate color from those DOF rays
            RenderState state {
                .scene = scene,
                .features = features,
                .bvh = bvh,
                .sampler = Sampler()
            };

            glm::vec3 pixelColor = renderRays(state, dofRays, /* depth = */ 0);

            // Write final averaged color
            pixels[y * resolution.x + x] = pixelColor;
        }
    }
}



std::vector<Ray> generateDepthOfFieldRays(
    const Trackball& camera,
    const glm::vec2& pixelCenterNDC,
    const glm::ivec2& /* resolution */, // not used
    const Features& features,
    int sampleCount)
{
    std::vector<Ray> dofRays;
    dofRays.reserve(sampleCount);

    // Clear debug data for each new pixel
    debugData.lensPoints.clear();
    debugData.focalPoints.clear();
    debugData.imagePlanePoints.clear();

    // Get camera
    debugData.cameraPosition = camera.position();
    debugData.cameraForward = camera.forward();
    debugData.cameraUp = camera.up();

    // DOF settings
    float lensRadius = features.extra.dofSettings.lendRadius;
    float focalDistance = features.extra.dofSettings.focalDistance;
    float imagePlaneDist = features.extra.dofSettings.imagePlaneDistance;

    // Setup camera basis
    glm::vec3 forward = glm::normalize(camera.forward());
    glm::vec3 right = glm::normalize(glm::cross(forward, camera.up()));
    glm::vec3 up = glm::normalize(glm::cross(right, forward));

    // The "center" of the image plane is **in front** of the camera:
    glm::vec3 imagePlaneCenter = camera.position() + forward * imagePlaneDist;

    for (int s = 0; s < sampleCount; s++) {
        // 1) Sample random point on lens disk
        glm::vec2 lensSample = sampleDisk(lensRadius);
        glm::vec3 lensPoint = camera.position()
            + lensSample.x * right
            + lensSample.y * up;

        // 2) The point on the image plane for this pixel
        float ndcX = pixelCenterNDC.x; // [-1..1]
        float ndcY = pixelCenterNDC.y; // [-1..1]

        glm::vec3 imagePlanePoint = imagePlaneCenter
            + ndcX * right * imagePlaneDist
            + ndcY * up * imagePlaneDist;

        // 3) Direction from lens to that image-plane point
        glm::vec3 lensToImageDir = glm::normalize(imagePlanePoint - lensPoint);

        // 4) Extend that direction out by focalDistance to find the focal point
        glm::vec3 focalPoint = lensPoint + lensToImageDir * focalDistance;

        // 5) The final ray is lens -> focalPoint
        glm::vec3 rayDir = glm::normalize(focalPoint - lensPoint);
        dofRays.emplace_back(Ray(lensPoint, rayDir));

        // Save debug data
        debugData.lensPoints.push_back(lensPoint);
        debugData.focalPoints.push_back(focalPoint);
        debugData.imagePlanePoints.push_back(imagePlanePoint);
    }

    return dofRays;
}





// TODO; Extra feature
// Given the same input as for `renderImage()`, instead render an image with your own implementation
// of motion blur. Here, you integrate over a time domain, and not just the pixel's image domain,
// to give objects the appearance of "fast movement".
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderImageWithMotionBlur(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    if (!features.extra.enableMotionBlur) {
        return;
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution(0.0, 1.0);

    const BSpline& bSpline = scene.bSpline;
    Trackball interpolatedCamera { camera };

    for (int y = 0; y < screen.resolution().y; y++) {
        for (int x = 0; x != screen.resolution().x; x++) {
            glm::vec3 avgCol { 0 };
            for (int i = 0; i < bSpline.rayNumber; i++) {
                float t = distribution(generator);
                SplinePoint point = interpolateBSpline(bSpline, t);

                float yaw = point.rotation.x;
                float pitch = point.rotation.y;

                glm::vec3 direction = glm::vec3(-pitch, -yaw + glm::pi<float>() / 2, 0); 

                interpolatedCamera.setCamera(point.position, direction, 0.0f);

                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };

                auto rays = generatePixelRays(state, interpolatedCamera, { x, y }, screen.resolution());
                auto L = renderRays(state, rays);
                avgCol += L;
            }
            avgCol /= float(bSpline.rayNumber);
            screen.setPixel(x, y, avgCol);
        }
    }
}


// Helper for postprocessImageWithBloom
unsigned long long factorial(int num) {
    unsigned long long result = 1;
    for (int i = 1; i <= num; i++) {
        result *= i;
    }
    return result;
}


// TODO; Extra feature
// Given a rendered image, compute and apply a bloom post-processing effect to increase bright areas.
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void postprocessImageWithBloom(const Scene& scene, const Features& features, const Trackball& camera, Screen& image)
{
    if (!features.extra.enableBloomEffect) {
        return;
    }


    auto& pixels = image.pixels();
    glm::ivec2 resoultion = image.resolution();
    std::vector<glm::vec3> brightMask(pixels.size());

    float minThreshold = features.extra.bloom.minThreshold; 
    float maxThreshold = features.extra.bloom.maxThreshold;

    // Mapping options
    const char* mappingOptions[] = {
        "Binary",
        "Linear",
        "Exponential",
        "Logartihmich",
        "Sigmoid",
        "Piecewise"
    };

    MappingType selectedMapping = features.extra.bloom.mappingType;


    // Piecewise Mapping
    for (int i = 0; i < pixels.size(); i++) {
        float luminance = glm::dot(pixels[i], glm::vec3(0.2126f, 0.7152f, 0.0722f));
        float brightness = 0.0f;
        float normalized = 0.0f;      

        // Select Mapping
        switch (selectedMapping) {

        case MappingType::Binary: {
            if (luminance > minThreshold) {
                brightness = 1.0f;
            } else {
                brightness = 0.0f;
            }
        } break;

        case MappingType::Linear: {
            brightness = glm::clamp((luminance - minThreshold) / (maxThreshold - minThreshold), 0.0f, 1.0f);
        } break;

        case MappingType::Exponential: {
            normalized = glm::clamp((luminance - minThreshold) / (maxThreshold - minThreshold), 0.0f, 1.0f);
            brightness = std::pow(normalized, 2.0f);
        } break;

        case MappingType::Logarithmic: {
            normalized = glm::clamp((luminance - minThreshold) / (maxThreshold - minThreshold), 0.0f, 1.0f);
            brightness = std::log(1.0f + 10.0f * normalized);
        } break;

        case MappingType::Sigmoid: {
            float midPoint = (minThreshold + maxThreshold) / 2.0f;
            brightness = 1.0f / (1.0f + std::exp(-10.0f * (luminance - midPoint)));
        } break;

        case MappingType::Piecewise: {
            if (luminance < minThreshold) {
                brightness = 0.0f;
            } else if (luminance < ((minThreshold + maxThreshold) / 2.0f)) // mid-brightness
            {
                brightness = glm::clamp((luminance - minThreshold) / ((minThreshold + maxThreshold) / 2.0f - minThreshold), 0.0f, 1.0f);
            } else // high brightness
            {
                normalized = glm::clamp((luminance - minThreshold) / (maxThreshold - minThreshold), 0.0f, 1.0f);
                brightness = std::pow(normalized, 2.0f);
            }
        } break;

        default:
            brightness = 0.0f;
            break;
        }

        brightMask[i] = pixels[i] * brightness;
    }

    // 1D Gaussian kernel using binomial
    const int kernelSize = features.extra.bloom.kernelSize; // TODO implement interface
    std::vector<float> kernel(kernelSize);
    int n = kernelSize - 1;
    float sum = 0.0f;

    // Binomial Coefficents
    for (int k = 0; k <= n; k++) 
    {
        kernel[k] = factorial(n) / (factorial(k) * factorial(n - k));
        sum += kernel[k];
    }

    for (auto& value : kernel)
    {
        value /= sum; // Normalize
    }

    // Horizontal 
    std::vector<glm::vec3> tempMask = brightMask;
    for (int y = 0; y < resoultion.y; y++) {
        for (int x = 0; x < resoultion.x; x++) {
            glm::vec3 blurredPixel(0.0f);

            for (int k = -kernelSize / 2; k <= kernelSize / 2; k++) {
                int neighbourIndex = glm::clamp(x + k, 0, resoultion.x - 1);
                blurredPixel += brightMask[y * resoultion.x + neighbourIndex] * kernel[k + kernelSize / 2];
            }
            tempMask[y * resoultion.x + x] = blurredPixel;
        }
    }

    // Vertical
    for (int x = 0; x < resoultion.x; x++) {
        for (int y = 0; y < resoultion.y; y++) {
            glm::vec3 blurredPixel(0.0f);
            for (int k = -kernelSize / 2; k <= kernelSize / 2; k++) {
                int neighbourIndex = glm::clamp(y + k, 0, resoultion.y - 1);
                blurredPixel += tempMask[neighbourIndex * resoultion.x + x] * kernel[k + kernelSize / 2];
            }
            brightMask[y * resoultion.x + x] = blurredPixel;
        }
    }

    // Combine with original image
    for (int i = 0; i < pixels.size(); i++)
    {
        pixels[i] += brightMask[i]; // Add bloom
    }

    // Shows only the Bloom contribution
    if (features.extra.bloom.onlyBloom)
    {
        for (int i = 0; i < pixels.size(); i++)
        {
            pixels[i] = brightMask[i];

            if (features.extra.bloom.onlyBloomGrayscale) 
            {
                float luminance = glm::dot(pixels[i], glm::vec3(0.2126f, 0.7152f, 0.0722f));
                pixels[i] = glm::vec3(luminance);
            }
        }
    }




}

void drawGlossyRayDebug(std::vector<glm::vec2> offsetsInCircle);
void drawOrthogonal(Ray r, glm::vec3 baseVector1, glm::vec3 baseVector2);

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) and an intersection, evaluates the contribution of a set of
// glossy reflective rays, recursively evaluating renderRay(..., depth + 1) along each ray, and adding the
// results times material.ks to the current intersection's hit color.
// - state;    the active scene, feature config, bvh, and sampler
// - ray;      camera ray
// - hitInfo;  intersection object
// - hitColor; current color at the current intersection, which this function modifies
// - rayDepth; current recursive ray depth
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
void renderRayGlossyComponent(RenderState& state, Ray ray, const HitInfo& hitInfo, glm::vec3& hitColor, int rayDepth)
{
    // Generate an initial specular ray, and base secondary glossies on this ray
    // auto numSamples = state.features.extra.numGlossySamples;
    // ...

    int numSamples = state.features.extra.numGlossySamples;

    Ray r = generateReflectionRay(ray, hitInfo);

    // compute the two base vectors of the disk
    glm::vec3 baseVector1 = glm::vec3(0, 0, 0);
    if (r.direction.x != 0) {
        baseVector1.y = 1;
        baseVector1.z = 1;
        baseVector1.x = (- r.direction.y - r.direction.z) / r.direction.x;
    } else if (r.direction.y != 0) {
        baseVector1.x = 1;
        baseVector1.z = 1;
        baseVector1.y = (- r.direction.x - r.direction.z) / r.direction.y;
    } else {
        baseVector1.x = 1;
        baseVector1.y = 1;
        baseVector1.z = (- r.direction.x - r.direction.y) / r.direction.z;
    }
    baseVector1 = normalize(baseVector1);
    glm::vec3 baseVector2 = normalize(cross(r.direction, baseVector1));

    float radius = 0.5f;
    // store the offsets so that we can draw the samples from the circle
    std::vector<glm::vec2> offsetsInCircle;

    // sample the rays in the disk
    for (int i = 0; i < numSamples; i++) {
        float angle = state.sampler.next_1d() * 2.0f * glm::pi<float>();
        Ray glossyRay = r;
        float radiusOffset = state.sampler.next_1d();
        float offset1 = radiusOffset *  sin(angle);
        float offset2 = radiusOffset * cos(angle);
        // add the offsets(y is inverted as the UI window postions (0, 0) at top left corner)
        offsetsInCircle.push_back(glm::vec2(offset1, -offset2));
        glossyRay.direction = r.direction + offset1 * radius * baseVector1 + offset2 * radius * baseVector2;

        // prevent rays from hitting the current surface again
        if (dot(glossyRay.direction, hitInfo.normal) <= 0) {
            continue;
        }

        // add the color of the glossy rays
        hitColor += renderRay(state, glossyRay, rayDepth + 1) * hitInfo.material.ks * glm::vec3(1.0f / numSamples);
    }

    // draw the distribution of the samples(only for the first reflection)
    if (state.features.enableDebugDraw && rayDepth == 0) {
        drawGlossyRayDebug(offsetsInCircle);
        drawOrthogonal(r, baseVector1, baseVector2);
    }
}

void drawOrthogonal(Ray r, glm::vec3 baseVector1, glm::vec3 baseVector2) {
    r.t = 1.0f;
    drawRay(r, glm::vec3(1, 0, 1));

    Ray baseRay1 = r;
    baseRay1.direction = baseVector1;
    drawRay(baseRay1, glm::vec3(1, 0, 1));

    Ray baseRay2 = r;
    baseRay2.direction = baseVector2;
    drawRay(baseRay2, glm::vec3(1, 0, 1));
}

void drawGlossyRayDebug(std::vector<glm::vec2> offsetsInCircle) {
    ImGui::Begin("Glossy Reflection Visual Debug");

    ImVec2 imagePos = ImGui::GetCursorScreenPos();
    ImDrawList* drawList = ImGui::GetWindowDrawList();

    drawList->AddCircle(ImVec2(imagePos.x + 150, imagePos.y + 150), 150.0f, IM_COL32(255, 0, 0, 255));
    drawList->AddCircleFilled(ImVec2(imagePos.x + 150, imagePos.y + 150), 5.0f, IM_COL32(255, 0, 0, 255));

    // draw the random samples and calculate the spread;
    float spreadSumX = 0.0f;
    float spreadSumY = 0.0f;
    for (int i = 0; i < offsetsInCircle.size(); i++) {
        drawList->AddCircleFilled(ImVec2(imagePos.x + 150 + offsetsInCircle[i].x * 150,
         imagePos.y + 150 + offsetsInCircle[i].y * 150),
         3.0f, IM_COL32(255, 255, 255, 255));
         spreadSumX += offsetsInCircle[i].x;
         spreadSumY += -offsetsInCircle[i].y;
    }
    ImGui::Dummy(ImVec2(0.0f, 300.0f));
    // evaluation how spread are the samples
    ImGui::Text("Deviation in the x direction: %.3f", spreadSumX / offsetsInCircle.size());
    ImGui::Text("Deviation in the y direction: %.3f", spreadSumY / offsetsInCircle.size());
    ImGui::End();
}

// TODO; Extra feature
// Given a camera ray (or reflected camera ray) that does not intersect the scene, evaluates the contribution
// along the ray, originating from an environment map. You will have to add support for environment textures
// to the Scene object, and provide a scene with the right data to supply this.
// - state; the active scene, feature config, bvh, and sampler
// - ray;   ray object
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
glm::vec3 sampleEnvironmentMap(RenderState& state, Ray ray)
{
    if (state.features.extra.enableEnvironmentMap) {
        // Part of your implementation should go here

        float mapSize = state.scene.mapSize;
        std::vector<Mesh> environmentMeshes = state.scene.environmentMap;

        // check which environment map we need
        if (state.scene.isSpherePlanet) {
            environmentMeshes = state.scene.environmentMapSpherePlanet;
        }
        if (state.scene.isCubePlanet) {
            environmentMeshes = state.scene.environmentMapCubePlanet;
        }

        // regulate the size of the environment map
        for (int i = 0; i < environmentMeshes.size(); i++) {
            for (int j = 0; j < environmentMeshes[i].vertices.size(); j++) {
                environmentMeshes[i].vertices[j].position *= mapSize;
            }
        }

        if (state.features.enableDebugDraw) {
            for (const auto& mesh : environmentMeshes) {
                drawMesh(mesh);
            }
        }

        HitInfo hitInfo;

        bool hit = false;

        // check if we have maximum size - simulate infinite;y far away environment map
        if (mapSize == 64) {
            ray.origin = glm::vec3(0.0f);
        }

        glm::vec3 hitPosition;

        // compute the texture from the environment map
        for (int i = 0; i < environmentMeshes.size(); i++) {
            Mesh mesh = environmentMeshes[i];
            for (int j = 0; j < mesh.triangles.size(); j++) {
                glm::uvec3 triangle = mesh.triangles[j];

                ray.t = FLT_MAX;
                if (!intersectRayWithTriangle(mesh.vertices[triangle.x].position, mesh.vertices[triangle.y].position,
                 mesh.vertices[triangle.z].position, ray, hitInfo)) {
                    continue;
                }

                hit = true;

                hitInfo.material = mesh.material;
                hitInfo.barycentricCoord = computeBarycentricCoord(mesh.vertices[triangle.x].position, mesh.vertices[triangle.y].position,
                 mesh.vertices[triangle.z].position, ray.origin + ray.direction * ray.t);
                hitInfo.texCoord = interpolateTexCoord(mesh.vertices[triangle.x].texCoord, mesh.vertices[triangle.y].texCoord,
                 mesh.vertices[triangle.z].texCoord, hitInfo.barycentricCoord);

                hitPosition = hitInfo.barycentricCoord.x * mesh.vertices[triangle.x].position + 
                hitInfo.barycentricCoord.y * mesh.vertices[triangle.y].position + 
                hitInfo.barycentricCoord.z * mesh.vertices[triangle.z].position;
                // glm::vec3 interpolatedNormal = interpolateNormal(
                //         mesh.vertices[triangle.x].normal,
                //         mesh.vertices[triangle.y].normal,
                //         mesh.vertices[triangle.z].normal,
                //         hitInfo.barycentricCoord);
                // hitInfo.normal = glm::normalize(interpolatedNormal);
                break;
            }
        }

        if (!hit) {
            return glm::vec3(0);
        }

        glm::vec2 hitPos = glm::vec2(hitInfo.texCoord.y, hitInfo.texCoord.x);

        glm::vec3 result = sampleTextureBilinear(*hitInfo.material.kdTexture, hitPos);

        if (state.features.enableDebugDraw) {
            drawRay(ray, result);
            drawSphere(hitPosition, 0.5f, result);
        }

        return result;
    } else {
        return glm::vec3(0.f);
    }
}
