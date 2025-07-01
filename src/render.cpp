#include "render.h"
#include "bvh_interface.h"
#include "draw.h"
#include "extra.h"
#include "light.h"
#include "recursive.h"
#include "sampler.h"
#include "screen.h"
#include "shading.h"
#include <framework/trackball.h>

#include <GLFW/glfw3.h>
#include <iostream>

#ifdef NDEBUG
#include <omp.h>
#endif


// This function is provided as-is. You do not have to implement it.
// Given relevant objects (scene, bvh, camera, etc) and an output screen, multithreaded fills
// each of the pixels using one of the below `renderPixel*()` functions, dependent on scene
// configuration. By default, `renderPixelNaive()` is called.
void renderImage(const Scene& scene, const BVHInterface& bvh, const Features& features, const Trackball& camera, Screen& screen)
{
    // Either directly render the image, or pass through to extra.h methods
    if (features.extra.enableDepthOfField && !features.enableDebugDraw) {
        renderImageWithDepthOfField(scene, bvh, features, camera, screen);
    } else if (features.extra.enableMotionBlur) {
        renderImageWithMotionBlur(scene, bvh, features, camera, screen);
    } else {
#ifdef NDEBUG // Enable multi threading in Release mode
#pragma omp parallel for schedule(guided)
#endif
        for (int y = 0; y < screen.resolution().y; y++) {
            for (int x = 0; x != screen.resolution().x; x++) {
                // Assemble useful objects on a per-pixel basis; e.g. a per-thread sampler
                // Note; we seed the sampler for consistenct behavior across frames
                RenderState state = {
                    .scene = scene,
                    .features = features,
                    .bvh = bvh,
                    .sampler = { static_cast<uint32_t>(screen.resolution().y * x + y) }
                };
                std::vector<Ray> rays;
               
                // Fall back to default ray generation
                rays = generatePixelRays(state, camera, { x, y }, screen.resolution());
            
                auto L = renderRays(state, rays);
                screen.setPixel(x, y, L);
            }
        }
    }

    // Pass through to extra.h for post processing
    if (features.extra.enableBloomEffect) {
        postprocessImageWithBloom(scene, features, camera, screen);
    }
}

// This function is provided as-is. You do not have to implement it.
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples for this pixel.
// This method forwards to `generatePixelRaysMultisampled` and `generatePixelRaysStratified` when necessary.
std::vector<Ray> generatePixelRays(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    std::vector<Ray> rays;

    // Check if depth-of-field is enabled
    if (state.features.extra.enableDepthOfField) {
        // Pixel center in NDC
        glm::vec2 pixelCenter = glm::vec2(
            (pixel.x + 0.5f) / screenResolution.x * 2.0f - 1.0f,
            1.0f - (pixel.y + 0.5f) / screenResolution.y * 2.0f);

        // Call generateDepthOfFieldRays to create rays for this pixel
        rays = generateDepthOfFieldRays(
            camera,
            pixelCenter,
            screenResolution,
            state.features,
            state.features.numPixelSamples);

        std::cout << "Rays are generated" << std::endl;
    } else {
        // Fall back to default ray generation
        if (state.features.numPixelSamples > 1) {
            if (state.features.enableJitteredSampling) {
                return generatePixelRaysStratified(state, camera, pixel, screenResolution);
            } else {
                return generatePixelRaysUniform(state, camera, pixel, screenResolution);
            }
        } else {
            // Generate single camera ray placed at the pixel's center
            glm::vec2 position = (glm::vec2(pixel) + 0.5f) / glm::vec2(screenResolution) * 2.f - 1.f;
            return { camera.generateRay(position) };
        }
    }

    return rays;
}

void renderPixelGridWindow(std::vector<glm::vec2> rayPositions, int numSamples);
void drawPixelGrid(std::vector<glm::vec2> rayPositions, int numSamples);

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// uniformly throughout this pixel.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is unit-tested, so do not change the function signature.
std::vector<Ray> generatePixelRaysUniform(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples camera rays uniformly distributed across the pixel. Use
    // Hint; use `state.sampler.next*d()` to generate random samples in [0, 1).
    auto numSamples = state.features.numPixelSamples;
    std::vector<Ray> rays;
    // ...

    std::vector<glm::vec2> rayPositions;

    for (int i = 0; i < numSamples; i++) {
        glm::vec2 randomPos = state.sampler.next_2d();
        rayPositions.push_back(glm::vec2(randomPos.x * 2 - 1.0f, randomPos.y * 2 - 1.0f));
        glm::vec2 rayPos = (glm::vec2(pixel) + randomPos) / glm::vec2(screenResolution) * 2.f - 1.f;
        rays.push_back(camera.generateRay(rayPos));
    }

    if (state.features.enableDebugDraw) {
        renderPixelGridWindow(rayPositions, std::round(std::sqrt(float(numSamples))));
    }

    return rays;
}

// TODO: standard feature
// Given a render state, camera, pixel position, and output resolution, generates a set of camera ray samples placed
// using jittered sampling throughout this pixel. Given NxN cells across one pixel, each ray sample is randomly
// placed somewhere within a cell.
// - state;            the active scene, feature config, bvh, and sampler
// - camera;           the camera object, used for ray generation
// - pixel;            x/y coordinates of the current pixel
// - screenResolution; x/y dimensions of the output image
// - return;           a vector of camera rays into the pixel
// This method is not unit-tested, but we do expect to find it **exactly here**, and we'd rather
// not go on a hunting expedition for your implementation, so please keep it here!
std::vector<Ray> generatePixelRaysStratified(RenderState& state, const Trackball& camera, glm::ivec2 pixel, glm::ivec2 screenResolution)
{
    // Generate numSamples * numSamples camera rays as jittered samples across the pixel.
    // Hint; use `state.sampler.next*d()` to generate random samples in [0, 1).
    auto numSamples = static_cast<uint32_t>(std::round(std::sqrt(float(state.features.numPixelSamples))));
    std::vector<Ray> rays;
    // ...

    std::vector<glm::vec2> rayPositions;

    for (int i = 0; i < numSamples; i++) {
        for (int j = 0; j < numSamples; j++) {
            float randomX = state.sampler.next_1d();
            float randomY = state.sampler.next_1d();

            rayPositions.push_back(glm::vec2(((float)i / numSamples + randomX / numSamples) * 2 - 1.0f,
             ((float)j / numSamples + randomY / numSamples) * 2 - 1.0f));

            glm::vec2 rayPos = glm::vec2(pixel.x + (float)i / numSamples + randomX / numSamples,
                pixel.y + (float)j / numSamples + randomY / numSamples) / 
                glm::vec2(screenResolution) * 2.f - 1.f;
            rays.push_back(camera.generateRay(rayPos));
        }
    }

   /*  if (state.features.enableDebugDraw) {
        renderPixelGridWindow(rayPositions, numSamples);
    }*/

    return rays;
}


void renderPixelGridWindow(std::vector<glm::vec2> rayPositions, int numSamples) {
    GLFWwindow* current = glfwGetCurrentContext();
    GLFWwindow* window = glfwCreateWindow(800, 600, "Pixel Grid", nullptr, nullptr);

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    while (!glfwWindowShouldClose(window)) {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        drawPixelGrid(rayPositions, numSamples);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwMakeContextCurrent(current);
}

void drawPixelGrid(std::vector<glm::vec2> rayPositions, int numSamples) {
    glColor3f(1.0f, 1.0f, 1.0f);

    float xPos = -1.0f + 2.0f / numSamples;
    for (int i = 0; i < numSamples; i++) {
        glBegin(GL_LINES);
        glVertex2f(xPos, -1.0f);
        glVertex2f(xPos, 1.0f);
        glEnd();

        xPos += 2.0f / numSamples;
    }

    float yPos = -1.0f + 2.0f / numSamples;
    for (int i = 0; i < numSamples; i++) {
        glBegin(GL_LINES);
        glVertex2f(-1.0f, yPos);
        glVertex2f(1.0f, yPos);
        glEnd();

        yPos += 2.0f / numSamples;
    }

    for (int i = 0; i < rayPositions.size(); i++) {
        glColor3f(1.0f, 0.0f, 0.0f);
        glPointSize(10.0f);
        glBegin(GL_POINTS);
        glVertex2f(rayPositions[i].x, rayPositions[i].y);
        glEnd();
    }
}