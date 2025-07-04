#include "bvh.h"
#include "config.h"
#include "draw.h"
#include "light.h"
#include "recursive.h"
#include "render.h"
#include "sampler.h"
#include "screen.h"
#include "splines.h"
#include "extra.h"

// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <fmt/chrono.h>
#include <fmt/core.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <imgui/imgui.h>
#include <nativefiledialog/nfd.h>
DISABLE_WARNINGS_POP()
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <framework/imgui_helper.h>
#include <framework/trackball.h>
#include <framework/variant_helper.h>
#include <framework/window.h>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <thread>
#include <variant>
#include <extra.h>
#include <extra.cpp>

// This is the main application. The code in here does not need to be modified.
enum class ViewMode {
    Rasterization = 0,
    RayTracing = 1
};

int debugBVHLeafId = 0;
uint32_t debugRaySeed = 4; // Chosen by fair dice roll

static void setOpenGLMatrices(const Trackball& camera);
static void drawLightsOpenGL(const Scene& scene, const Trackball& camera, int selectedLight);
static void drawSceneOpenGL(const Scene& scene);
bool sliderIntSquarePower(const char* label, int* v, int v_min, int v_max);

int main(int argc, char** argv)
{
    Config config = {};
    if (argc > 1) {
        config = readConfigFile(argv[1]);
    } else {
        // Add a default camera if no config file is given.
        config.cameras.emplace_back(CameraConfig {});
    }

    if (!config.cliRenderingEnabled) {
        Trackball::printHelp();
        std::cout << "\n Press the [R] key on your keyboard to create a ray towards the mouse cursor" << std::endl
                  << std::endl;

        Window window { "Final Project", config.windowSize, OpenGLVersion::GL2, true };
        Screen screen { config.windowSize, true };
        Trackball camera { &window, glm::radians(config.cameras[0].fieldOfView), config.cameras[0].distanceFromLookAt };
        camera.setCamera(config.cameras[0].lookAt, glm::radians(config.cameras[0].rotation), config.cameras[0].distanceFromLookAt);

        SceneType sceneType { SceneType::SingleTriangle };
        std::vector<Ray> debugRays;

        Scene scene = loadScenePrebuilt(sceneType, config.dataPath);
        BVH bvh(scene, config.features);

        int bvhDebugLevel = 0;
        int bvhDebugLeaf = 0;
        bool debugBVHLevel { false };
        bool debugBVHLeaf { false };
        ViewMode viewMode { ViewMode::Rasterization };

        window.registerKeyCallback([&](int key, int /* scancode */, int action, int /* mods */) {
            if (action == GLFW_PRESS) {
                switch (key) {
                case GLFW_KEY_R: {
                    // Shoot a ray. Produce a ray from camera to the far plane.
                    RenderState state = { .scene = scene, .features = config.features, .bvh = bvh, .sampler = { debugRaySeed } };
                    auto tmp = window.getNormalizedCursorPos();
                    auto pixel = glm::ivec2(tmp * glm::vec2(screen.resolution()));
                    debugRays = generatePixelRays(state, camera, pixel, screen.resolution());

                    /// genatrtdepth 
                } break;
                case GLFW_KEY_A: {
                    debugBVHLeafId++;
                } break;
                case GLFW_KEY_S: {
                    debugBVHLeafId = std::max(0, debugBVHLeafId - 1);

                } break;
                case GLFW_KEY_ESCAPE: {
                    window.close();
                } break;
                };
            }
        });

        const size_t points = 3;
        std::vector<glm::vec3> position(points);
        std::vector<glm::vec2> rotation(points);
        bool test = false;
  
        scene.bSpline = BSpline{position, rotation, false, Linear, 0, 1};




        int selectedLightIdx = scene.lights.empty() ? -1 : 0;
        while (!window.shouldClose()) {
            window.updateInput();

            // === Setup the UI ===
            ImGui::Begin("Final Project");
            {
                constexpr std::array items {
                    "SingleTriangle",
                    "Cube (segment light)",
                    "Cube (textured)",
                    "Cornell Box (with mirror)",
                    "Cornell Box (with transparency)",
                    "Cornell Box (parallelogram light and mirror)",
                    "Monkey",
                    "Teapot",
                    "Dragon",
                    /* "AABBs",*/ "Spheres", /*"Mixed",*/
                    "Custom",
                    "Environment map"
                };
                if (ImGui::Combo("Scenes", reinterpret_cast<int*>(&sceneType), items.data(), int(items.size()))) {
                    debugRays.clear();
                    scene = loadScenePrebuilt(sceneType, config.dataPath);
                    selectedLightIdx = scene.lights.empty() ? -1 : 0;
                    bvh = BVH(scene, config.features);

                    if (!debugRays.empty()) {
                        RenderState state = { .scene = scene, .features = config.features, .bvh = bvh, .sampler = { debugRaySeed } };
                        renderRays(state, debugRays);
                    }
                }
            }
            {
                constexpr std::array items { "Rasterization", "Ray Traced" };
                ImGui::Combo("View mode", reinterpret_cast<int*>(&viewMode), items.data(), int(items.size()));
            }

            ImGui::Separator();
            if (ImGui::CollapsingHeader("Features", ImGuiTreeNodeFlags_DefaultOpen)) {
                // Feature toggles
                ImGui::Checkbox("BVH", &config.features.enableAccelStructure);
                ImGui::Checkbox("Reflections", &config.features.enableReflections);
                ImGui::Checkbox("Normal interpolation", &config.features.enableNormalInterp);
                ImGui::Checkbox("Shading", &config.features.enableShading);
                if (config.features.enableShading) {
                    ImGui::Indent();
                    constexpr std::array items {
                        "Lambertian", "Phong", "Blinn-Phong", "Linear Gradient", "Linear Gradient Comparison"
                    };
                    if (ImGui::Combo("Shading model", reinterpret_cast<int*>(&config.features.shadingModel), items.data(), int(items.size()))) {
                        // Do nothing
                    }
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Shadows", &config.features.enableShadows);
                if (config.features.enableShadows) {
                    uint32_t minSamples = 1u, maxSamples = 64u;
                    ImGui::Indent();
                    ImGui::SliderScalar("Shadow samples", ImGuiDataType_U32, &config.features.numShadowSamples, &minSamples, &maxSamples);
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Textures", &config.features.enableTextureMapping);
                if (config.features.enableTextureMapping) {
                    ImGui::Indent();
                    ImGui::Checkbox("Texture filtering (bilinear)", &config.features.enableBilinearTextureFiltering);
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Transparency", &config.features.enableTransparency);

                ImGui::Separator(); // Pixel sampling does not have an accompanying toggle
                {
                    ImGui::Checkbox("Jittered sampling", &config.features.enableJitteredSampling);
                    uint32_t minSamples = 1u, maxSamples = 64u;
                    ImGui::SliderScalar("Pixel samples", ImGuiDataType_U32, &config.features.numPixelSamples, &minSamples, &maxSamples);
                }
            }

            if (ImGui::CollapsingHeader("Extra Features")) {
                ImGui::Checkbox("Bloom effect", &config.features.extra.enableBloomEffect);
                if (config.features.extra.enableBloomEffect) {
                    ImGui::Indent();
                    // Sliders for bloom settings
                    ImGui::SliderFloat("Min Threshold", &config.features.extra.bloom.minThreshold, 0.0f, 1.0f, "%.2f");
                    ImGui::SliderFloat("Max Threshold", &config.features.extra.bloom.maxThreshold, 0.0f, 1.0f, "%.2f");

                    // Slider for kernel size
                    if (ImGui::SliderInt("Kernel Size", &config.features.extra.bloom.kernelSize, 3, 25)) {
                        // Ensure kernel size is odd
                        if (config.features.extra.bloom.kernelSize % 2 == 0) {
                            config.features.extra.bloom.kernelSize++;
                        }
                    }

                     // Dropdown for mapping type
                    const char* mappingOptions[] = { "Binary", "Linear", "Exponential", "Logarithmic", "Sigmoid", "Piecewise" };
                    int currentMapping = static_cast<int>(config.features.extra.bloom.mappingType);
                    if (ImGui::Combo("Mapping Type", &currentMapping, mappingOptions, IM_ARRAYSIZE(mappingOptions))) {
                        config.features.extra.bloom.mappingType = MappingType(currentMapping);
                    }

                    // Checkbox for showing only the bloom contribution
                    ImGui::Checkbox("Show only bloom contribution", &config.features.extra.bloom.onlyBloom);

                     // Add a nested option for grayscale visualization
                    if (config.features.extra.bloom.onlyBloom) {
                        ImGui::Indent();
                        ImGui::Checkbox("Show bloom contribution as grayscale", &config.features.extra.bloom.onlyBloomGrayscale);
                        ImGui::Unindent();
                    }

                    ImGui::Unindent();
                }
                ImGui::Checkbox("Depth of field", &config.features.extra.enableDepthOfField);
                if (config.features.extra.enableDepthOfField) {
                    ImGui::Indent();
                    // Sliders for DOF settings
                    ImGui::SliderFloat("Lens Radius", &config.features.extra.dofSettings.lendRadius, 0.01f, 1.0f, "%.2f");
                    ImGui::SliderFloat("Focal Distance", &config.features.extra.dofSettings.focalDistance, 1.0f, 100.0f, "%.2f");
                    ImGui::SliderFloat("Image Plane Distance", &config.features.extra.dofSettings.imagePlaneDistance, 0.1f, 10.0f, "%.2f");
                    ImGui::Unindent();
                }


                ImGui::Checkbox("Motion blur", &config.features.extra.enableMotionBlur);
                if (config.features.extra.enableMotionBlur) {
                    ImGui::Indent();
                    // Add motion blur settings here, if necessary

                    constexpr std::array items {
                        "Linear",
                        "Quadratic",
                        "Cubic",
                    };

                    ImGui::Combo("Spline type", reinterpret_cast<int*>(&scene.bSpline.type), items.data(), int(items.size()));
                    ImGui::Checkbox("Connect endpoints", &scene.bSpline.repeatKnots);
                    ImGui::SliderFloat("View at t", &scene.bSpline.t, 0.0f, 1.0f);

                    
                    ImGui::Checkbox("Time point", &test);

                    if (test) {
                        SplinePoint point = interpolateBSpline(scene.bSpline, scene.bSpline.t);
                        glm::vec3 direction = glm::vec3(-point.rotation.y, -point.rotation.x + glm::pi<float>() / 2, 0); 
                        camera.setCamera(point.position, direction, 0.1f);
                    }
                    

                    size_t min_size = 2;
                    if (scene.bSpline.type == Quadratic)
                        min_size = 3;
                    else if (scene.bSpline.type == Cubic)
                        min_size = 4;

                    if (ImGui::Button("+")) {
                        scene.bSpline.position.push_back(glm::vec3 { 0 });
                        scene.bSpline.rotation.push_back(glm::vec2 { 0 });
                    }
                    if (ImGui::Button("-") && scene.bSpline.position.size() > min_size) {
                        scene.bSpline.position.pop_back();
                        scene.bSpline.rotation.pop_back();
                    }
                    
                    while (min_size > (int)scene.bSpline.position.size()) {
                        scene.bSpline.position.push_back(glm::vec3 { 0 });
                        scene.bSpline.rotation.push_back(glm::vec2 { 0 });
                    }

                   
                    for (int i = 0; i < scene.bSpline.position.size(); i++) {
                        std::string label = "Point " + std::to_string(i + 1);
                        ImGui::DragFloat3(label.c_str(), glm::value_ptr(scene.bSpline.position[i]), 0.01f, -3.0f, 3.0f, "%0.2f");

                        std::string label1 = "Yaw/Pitch" + std::to_string(i + 1);
                        ImGui::DragFloat2(label1.c_str(), glm::value_ptr(scene.bSpline.rotation[i]), 0.05f, -2 * glm::pi<float>(), 2.0 * glm::pi<float>(), "%0.2f");
                        //ImGui::InputFloat3("LookAt", glm::value_ptr(lookAt), "%0.2f", ImGuiInputTextFlags_ReadOnly);
                    }

                    ImGui::SliderInt("Rays repeated", &scene.bSpline.rayNumber, 1, 32);

                    ImGui::Unindent();
                }
                ImGui::Checkbox("Glossy reflections", &config.features.extra.enableGlossyReflection);
                if (config.features.extra.enableGlossyReflection) {
                    uint32_t minSamples = 1u, maxSamples = 64u;
                    ImGui::Indent();
                    ImGui::SliderScalar("Glossy samples", ImGuiDataType_U32, &config.features.extra.numGlossySamples, &minSamples, &maxSamples);
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Environment maps", &config.features.extra.enableEnvironmentMap);
                if (config.features.extra.enableEnvironmentMap) {
                    uint32_t minSamples = 1u, maxSamples = 64u;
                    ImGui::Indent();
                    static int selectedPlanetMap = 0;
                    if (ImGui::RadioButton("Normal map", selectedPlanetMap == 0)) {
                        selectedPlanetMap = 0;
                    }
                    if (ImGui::RadioButton("Sphere Planet Map", selectedPlanetMap == 1)) {
                        selectedPlanetMap = 1;
                    }
                        
                    if (ImGui::RadioButton("Cube Planet Map", selectedPlanetMap == 2)) {
                        selectedPlanetMap = 2;
                    }
                    scene.isSpherePlanet = (selectedPlanetMap == 1);
                    scene.isCubePlanet = (selectedPlanetMap == 2);
                    ImGui::SliderScalar("Environment size", ImGuiDataType_U32, &scene.mapSize, &minSamples, &maxSamples);
                    ImGui::Unindent();
                }
                ImGui::Checkbox("Texture filtering (mipmap)", &config.features.extra.enableMipmapTextureFiltering);
            }

            if (ImGui::TreeNode("Camera(read only)")) {
                auto lookAt = camera.lookAt();
                auto position = camera.position();
                auto rotation = glm::degrees(camera.rotationEulerAngles());
                auto distance = camera.distanceFromLookAt();
                ImGui::InputFloat3("Position", glm::value_ptr(position), "%0.2f", ImGuiInputTextFlags_ReadOnly);
                ImGui::InputFloat3("LookAt", glm::value_ptr(lookAt), "%0.2f", ImGuiInputTextFlags_ReadOnly);
                ImGui::InputFloat("Distance from look at", &distance, 0.1f, 0.1f, "%0.2f", ImGuiInputTextFlags_ReadOnly);
                ImGui::InputFloat3("Rotation", glm::value_ptr(rotation), "%0.2f", ImGuiInputTextFlags_ReadOnly);
                ImGui::TreePop();
            }

            ImGui::Spacing();
            ImGui::Separator();
            if (ImGui::Button("Render to file")) {
                // Show a file picker.
                nfdchar_t* pOutPath = nullptr;
                const nfdresult_t result = NFD_SaveDialog("bmp", nullptr, &pOutPath);
                if (result == NFD_OKAY) {
                    std::filesystem::path outPath { pOutPath };
                    free(pOutPath); // NFD is a C API so we have to manually free the memory it allocated.
                    outPath.replace_extension("bmp"); // Make sure that the file extension is *.bmp

                    // Perform a new render and measure the time it took to generate the image.
                    using clock = std::chrono::high_resolution_clock;
                    const auto start = clock::now();
                    renderImage(scene, bvh, config.features, camera, screen);
                    const auto end = clock::now();
                    std::cout << "Time to render image: " << std::chrono::duration<float, std::milli>(end - start).count() << " milliseconds" << std::endl;
                    // Store the new image.
                    screen.writeBitmapToFile(outPath);
                }
            }

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Text("Debugging");
            if (viewMode == ViewMode::Rasterization) {
                ImGui::Checkbox("Enable Debug Draw", &config.features.enableDebugDraw);
                ImGui::Checkbox("Draw BVH Level", &debugBVHLevel);
                if (debugBVHLevel)
                    ImGui::SliderInt("BVH Level", &bvhDebugLevel, 0, bvh.numLevels() - 1);
                ImGui::Checkbox("Draw BVH Leaf", &debugBVHLeaf);
                if (debugBVHLeaf)
                    ImGui::SliderInt("BVH Leaf", &bvhDebugLeaf, 1, bvh.numLeaves());
            }

            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Text("Lights");
            {
                std::vector<std::string> options;
                options.push_back("None");
                for (size_t i = 0; i < scene.lights.size(); i++) {
                    options.push_back("Light " + std::to_string(i));
                }
                std::vector<const char*> optionsPointers;
                std::transform(std::begin(options), std::end(options), std::back_inserter(optionsPointers), [](const auto& str) { return str.c_str(); });

                // Offset such that selectedLightIdx=-1 becomes item 0 (None).
                ++selectedLightIdx;
                ImGui::Combo("Selected light", &selectedLightIdx, optionsPointers.data(), static_cast<int>(optionsPointers.size()));
                --selectedLightIdx;

                if (selectedLightIdx >= 0) {
                    setOpenGLMatrices(camera);
                    std::visit(
                        make_visitor(
                            [&](PointLight& light) {
                                showImGuizmoTranslation(window, camera, light.position); // 3D controls to translate light source.
                                ImGui::DragFloat3("Light position", glm::value_ptr(light.position), 0.01f, -3.0f, 3.0f);
                                ImGui::ColorEdit3("Light color", glm::value_ptr(light.color), ImGuiColorEditFlags_Float | ImGuiColorEditFlags_HDR);
                            },
                            [&](SegmentLight& light) {
                                static int selectedEndpoint = 0;
                                // 3D controls to translate light source.
                                if (selectedEndpoint == 0)
                                    showImGuizmoTranslation(window, camera, light.endpoint0);
                                else
                                    showImGuizmoTranslation(window, camera, light.endpoint1);

                                const std::array<const char*, 2> endpointOptions { "Endpoint 0", "Endpoint 1" };
                                ImGui::Combo("Selected endpoint", &selectedEndpoint, endpointOptions.data(), int(endpointOptions.size()));
                                ImGui::DragFloat3("Endpoint 0", glm::value_ptr(light.endpoint0), 0.01f, -3.0f, 3.0f);
                                ImGui::DragFloat3("Endpoint 1", glm::value_ptr(light.endpoint1), 0.01f, -3.0f, 3.0f);
                                ImGui::ColorEdit3("Color 0", glm::value_ptr(light.color0));
                                ImGui::ColorEdit3("Color 1", glm::value_ptr(light.color1));
                            },
                            [&](ParallelogramLight& light) {
                                glm::vec3 vertex1 = light.v0 + light.edge01;
                                glm::vec3 vertex2 = light.v0 + light.edge02;

                                static int selectedVertex = 0;
                                // 3D controls to translate light source.
                                if (selectedVertex == 0)
                                    showImGuizmoTranslation(window, camera, light.v0);
                                else if (selectedVertex == 1)
                                    showImGuizmoTranslation(window, camera, vertex1);
                                else
                                    showImGuizmoTranslation(window, camera, vertex2);

                                const std::array<const char*, 3> vertexOptions { "Vertex 0", "Vertex 1", "Vertex 2" };
                                ImGui::Combo("Selected vertex", &selectedVertex, vertexOptions.data(), int(vertexOptions.size()));
                                ImGui::DragFloat3("Vertex 0", glm::value_ptr(light.v0), 0.01f, -3.0f, 3.0f);
                                ImGui::DragFloat3("Vertex 1", glm::value_ptr(vertex1), 0.01f, -3.0f, 3.0f);
                                light.edge01 = vertex1 - light.v0;
                                ImGui::DragFloat3("Vertex 2", glm::value_ptr(vertex2), 0.01f, -3.0f, 3.0f);
                                light.edge02 = vertex2 - light.v0;

                                ImGui::ColorEdit3("Color 0", glm::value_ptr(light.color0));
                                ImGui::ColorEdit3("Color 1", glm::value_ptr(light.color1));
                                ImGui::ColorEdit3("Color 2", glm::value_ptr(light.color2));
                                ImGui::ColorEdit3("Color 3", glm::value_ptr(light.color3));
                            },
                            [](auto) { /* any other type of light */ }),
                        scene.lights[size_t(selectedLightIdx)]);
                }
            }

            if (ImGui::Button("Add point light")) {
                selectedLightIdx = int(scene.lights.size());
                scene.lights.emplace_back(PointLight { .position = glm::vec3(0.0f), .color = glm::vec3(1.0f) });
            }
            if (ImGui::Button("Add segment light")) {
                selectedLightIdx = int(scene.lights.size());
                scene.lights.emplace_back(SegmentLight { .endpoint0 = glm::vec3(0.0f), .endpoint1 = glm::vec3(1.0f), .color0 = glm::vec3(1, 0, 0), .color1 = glm::vec3(0, 0, 1) });
            }
            if (ImGui::Button("Add parallelogram light")) {
                selectedLightIdx = int(scene.lights.size());
                scene.lights.emplace_back(ParallelogramLight {
                    .v0 = glm::vec3(0.0f),
                    .edge01 = glm::vec3(1, 0, 0),
                    .edge02 = glm::vec3(0, 1, 0),
                    .color0 = glm::vec3(1, 0, 0), // red
                    .color1 = glm::vec3(0, 1, 0), // green
                    .color2 = glm::vec3(0, 0, 1), // blue
                    .color3 = glm::vec3(1, 1, 1) // white
                });
            }
            if (selectedLightIdx >= 0 && ImGui::Button("Remove selected light")) {
                scene.lights.erase(std::begin(scene.lights) + selectedLightIdx);
                selectedLightIdx = -1;
            }

            // Clear screen.
            glViewport(0, 0, window.getFrameBufferSize().x, window.getFrameBufferSize().y);
            glClearDepth(1.0);
            glClearColor(0.0, 0.0, 0.0, 0.0);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            setOpenGLMatrices(camera);

            // Draw either using OpenGL (rasterization) or the ray tracing function.
            switch (viewMode) {
            case ViewMode::Rasterization: {
                glPushAttrib(GL_ALL_ATTRIB_BITS);
                if (debugBVHLeaf) {
                    glEnable(GL_POLYGON_OFFSET_FILL);
                    // To ensure that debug draw is always visible, adjust the scale used to calculate the depth value.
                    glPolygonOffset(float(1.4), 1.0);
                    drawSceneOpenGL(scene);
                    glDisable(GL_POLYGON_OFFSET_FILL);
                } else {
                    drawSceneOpenGL(scene);
                }

      

                {
                    if (!debugRays.empty()) {
                        // Call renderRay for the debug ray. Ignore the result but tell the function that it should
                        // draw the rays instead.
                        glDisable(GL_LIGHTING);
                        glDepthFunc(GL_LEQUAL);
                        RenderState state = { .scene = scene, .features = config.features, .bvh = bvh, .sampler = { debugRaySeed } };
                        (void)renderRays(state, debugRays);
                    }
                }
                glPopAttrib();

                drawLightsOpenGL(scene, camera, selectedLightIdx);

                if (debugBVHLevel || debugBVHLeaf) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    setOpenGLMatrices(camera);
                    glDisable(GL_LIGHTING);
                    glEnable(GL_DEPTH_TEST);

                    // Enable alpha blending. More info at:
                    // https://learnopengl.com/Advanced-OpenGL/Blending
                    glEnable(GL_BLEND);
                    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                    if (debugBVHLevel)
                        bvh.debugDrawLevel(bvhDebugLevel);
                    if (debugBVHLeaf)
                        bvh.debugDrawLeaf(bvhDebugLeaf);
                    glPopAttrib();
                }

                if (config.features.extra.enableMotionBlur) {
                    for (size_t i = 0; i < scene.bSpline.position.size(); i++) {
                        drawSphere(scene.bSpline.position[i], 0.01f);
                        drayRayYawPitch(scene.bSpline.position[i], scene.bSpline.rotation[i]);
                    }
                    for (size_t i = 0; i < scene.bSpline.position.size() - 1 && scene.bSpline.position.size() > 0; i++) {
                        drawLine(scene.bSpline.position[i], scene.bSpline.position[i + 1]);
                    }
                    drawSpline(scene.bSpline);
                }
            } break;
            case ViewMode::RayTracing: {
                config.features.enableDebugDraw = false;
                screen.clear(glm::vec3(0.0f));

                using clock = std::chrono::high_resolution_clock;
                const auto start = clock::now();
                renderImage(scene, bvh, config.features, camera, screen);
                const auto end = clock::now();
                const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                fmt::print("Rendering took {} ms.\n", duration);
                screen.setPixel(0, 0, glm::vec3(1.0f));
                screen.draw(); // Takes the image generated using ray tracing and outputs it to the screen using OpenGL.
            } break;
            default:
                break;
            }

            ImGui::End();
            window.swapBuffers();
        }
    } else {
        // Command-line rendering.
        std::cout << config;
        // NOTE(Yang): Trackball is highly coupled with the window,
        // so we need to create a dummy window here but not show it.
        // In this case, GLEW will not be initialized, OpenGL functions
        // will not be loaded, and the window will not be shown.
        // All debug draw calls will be disabled.
        Window window { "Final Project", config.windowSize, OpenGLVersion::GL2, false };
        // Load scene.
        Scene scene;
        std::string sceneName;
        std::visit(make_visitor(
                       [&](const std::filesystem::path& path) {
                           scene = loadSceneFromFile(path, config.lights);
                           sceneName = path.stem().string();
                       },
                       [&](const SceneType& type) {
                           scene = loadScenePrebuilt(type, config.dataPath);
                           sceneName = serialize(type);
                       }),
            config.scene);

        BVH bvh(scene, config.features);

        using clock = std::chrono::high_resolution_clock;
        // Create output directory if it does not exist.
        if (!std::filesystem::exists(config.outputDir)) {
            std::filesystem::create_directories(config.outputDir);
        }
        const auto start = clock::now();
        std::string start_time_string = fmt::format("{:%Y-%m-%d_%H-%M-%S}", fmt::localtime(std::time(nullptr)));

        for (std::size_t i = 0; i < config.cameras.size(); ++i) {
            const auto& cameraConfig = config.cameras[i];
            Screen screen { config.windowSize, false };
            screen.clear(glm::vec3(0.0f));
            Trackball camera { &window, glm::radians(cameraConfig.fieldOfView), cameraConfig.distanceFromLookAt };
            camera.setCamera(cameraConfig.lookAt, glm::radians(cameraConfig.rotation), cameraConfig.distanceFromLookAt);
            renderImage(scene, bvh, config.features, camera, screen);
            const auto filename_base = fmt::format("{}_{}_cam_{}", sceneName, start_time_string, i);
            const auto filepath = config.outputDir / (filename_base + ".bmp");
            fmt::print("Image {} saved to {}\n", i, filepath.string());
            screen.writeBitmapToFile(filepath);
        }
        const auto end = clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        fmt::print("Rendering took {} ms, {} images rendered.\n", duration, config.cameras.size());
    }

    return 0;
}

static void setOpenGLMatrices(const Trackball& camera)
{
    // Load view matrix.
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    const glm::mat4 viewMatrix = camera.viewMatrix();
    glMultMatrixf(glm::value_ptr(viewMatrix));

    // Load projection matrix.
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    const glm::mat4 projectionMatrix = camera.projectionMatrix();
    glMultMatrixf(glm::value_ptr(projectionMatrix));
}

static void drawLightsOpenGL(const Scene& scene, const Trackball& camera, int /* selectedLight */)
{
    // Normals will be normalized in the graphics pipeline.
    glEnable(GL_NORMALIZE);
    // Activate rendering modes.
    glEnable(GL_DEPTH_TEST);
    // Draw front and back facing triangles filled.
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    // Interpolate vertex colors over the triangles.
    glShadeModel(GL_SMOOTH);

    glDisable(GL_LIGHTING);
    // Draw all non-selected lights.
    for (size_t i = 0; i < scene.lights.size(); i++) {
        std::visit(
            make_visitor(
                [](const PointLight& light) { drawSphere(light.position, 0.01f, light.color); },
                [](const SegmentLight& light) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glBegin(GL_LINES);
                    glColor3fv(glm::value_ptr(light.color0));
                    glVertex3fv(glm::value_ptr(light.endpoint0));
                    glColor3fv(glm::value_ptr(light.color1));
                    glVertex3fv(glm::value_ptr(light.endpoint1));
                    glEnd();
                    glPopAttrib();
                    drawSphere(light.endpoint0, 0.01f, light.color0);
                    drawSphere(light.endpoint1, 0.01f, light.color1);
                },
                [](const ParallelogramLight& light) {
                    glPushAttrib(GL_ALL_ATTRIB_BITS);
                    glBegin(GL_QUADS);
                    glColor3fv(glm::value_ptr(light.color0));
                    glVertex3fv(glm::value_ptr(light.v0));
                    glColor3fv(glm::value_ptr(light.color1));
                    glVertex3fv(glm::value_ptr(light.v0 + light.edge01));
                    glColor3fv(glm::value_ptr(light.color3));
                    glVertex3fv(glm::value_ptr(light.v0 + light.edge01 + light.edge02));
                    glColor3fv(glm::value_ptr(light.color2));
                    glVertex3fv(glm::value_ptr(light.v0 + light.edge02));
                    glEnd();
                    glPopAttrib();
                },
                [](auto) { /* any other type of light */ }),
            scene.lights[i]);
    }

    // Draw a colored sphere at the location at which the trackball is looking/rotating around.
    glDisable(GL_LIGHTING);
    drawSphere(camera.lookAt(), 0.01f, glm::vec3(0.2f, 0.2f, 1.0f));
}

void drawSceneOpenGL(const Scene& scene)
{
    // Activate the light in the legacy OpenGL mode.
    glEnable(GL_LIGHTING);

    // Tell OpenGL where the lights are (so it nows how to shade surfaces in the scene).
    // This is only used in the rasterization view. OpenGL only supports point lights so
    // we replace segment/parallelogram lights by point lights.
    int i = 0;
    const auto enableLight = [&](const glm::vec3& position, const glm::vec3 color) {
        const GLenum glLight = static_cast<GLenum>(GL_LIGHT0 + i);
        glEnable(glLight);
        const glm::vec4 position4 { position, 1 };
        glLightfv(glLight, GL_POSITION, glm::value_ptr(position4));
        const glm::vec4 color4 { glm::clamp(color, 0.0f, 1.0f), 1.0f };
        const glm::vec4 zero4 { 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(glLight, GL_AMBIENT, glm::value_ptr(zero4));
        glLightfv(glLight, GL_DIFFUSE, glm::value_ptr(color4));
        glLightfv(glLight, GL_SPECULAR, glm::value_ptr(zero4));
        // NOTE: quadratic attenuation doesn't work like you think it would in legacy OpenGL.
        // The distance is not in world space but in NDC space!
        glLightf(glLight, GL_CONSTANT_ATTENUATION, 1.0f);
        glLightf(glLight, GL_LINEAR_ATTENUATION, 0.0f);
        glLightf(glLight, GL_QUADRATIC_ATTENUATION, 0.0f);
        i++;
    };
    for (const auto& light : scene.lights) {
        std::visit(
            make_visitor(
                [&](const PointLight& pointLight) {
                    enableLight(pointLight.position, pointLight.color);
                },
                [&](const SegmentLight& segmentLight) {
                    // Approximate with two point lights: one at each endpoint.
                    enableLight(segmentLight.endpoint0, 0.5f * segmentLight.color0);
                    enableLight(segmentLight.endpoint1, 0.5f * segmentLight.color1);
                },
                [&](const ParallelogramLight& parallelogramLight) {
                    enableLight(parallelogramLight.v0, 0.25f * parallelogramLight.color0);
                    enableLight(parallelogramLight.v0 + parallelogramLight.edge01, 0.25f * parallelogramLight.color1);
                    enableLight(parallelogramLight.v0 + parallelogramLight.edge02, 0.25f * parallelogramLight.color2);
                    enableLight(parallelogramLight.v0 + parallelogramLight.edge01 + parallelogramLight.edge02, 0.25f * parallelogramLight.color3);
                },
                [](auto) { /* any other type of light */ }),
            light);
    }

    // Draw the scene and the ray (if any).
    drawScene(scene);
}

bool sliderIntSquarePower(const char* label, int* v, int v_min, int v_max)
{
    // Round to the nearest square power.
    const float v_rounded = std::round(std::sqrt(float(*v)));
    *v = static_cast<int>(v_rounded * v_rounded);
    return ImGui::SliderInt(label, v, v_min, v_max);
}