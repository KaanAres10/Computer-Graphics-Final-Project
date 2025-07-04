#pragma once
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <framework/mesh.h>
#include <framework/ray.h>
#include <optional>
#include <variant>
#include <vector>
#include "common.h"
#include "splines.h"

enum SceneType {
    SingleTriangle,
    Cube,
    CubeTextured,
    CornellBox,
    CornellBoxTransparency,
    CornellBoxParallelogramLight,
    Monkey,
    Teapot,
    Dragon,
    Spheres,
    Custom,
    ReflectiveSphere
};

struct Scene {
    using SceneLight = std::variant<PointLight, SegmentLight, ParallelogramLight>;

    SceneType type;
    std::vector<Mesh> meshes;
    std::vector<Sphere> spheres;
    std::vector<SceneLight> lights;

    // You can add your own objects (e.g. environment maps) here
    // ...
    std::vector<Mesh> environmentMap;
    std::vector<Mesh> environmentMapSpherePlanet;
    std::vector<Mesh> environmentMapCubePlanet;
    bool isSpherePlanet = false;
    bool isCubePlanet = false;
    int mapSize = 1;
    BSpline bSpline;
};

// Load a prebuilt scene.
Scene loadScenePrebuilt(SceneType type, const std::filesystem::path& dataDir);

// Load a scene from a file.
Scene loadSceneFromFile(const std::filesystem::path& path, const std::vector<std::variant<PointLight, SegmentLight, ParallelogramLight>>& lights);
