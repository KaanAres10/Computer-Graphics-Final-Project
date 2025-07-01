#include "texture.h"
#include "render.h"
#include <framework/image.h>
#include <fmt/core.h>
#include <glm/gtx/string_cast.hpp>

#include <imgui/imgui.h>
#include <GLFW/glfw3.h>
#include "scene.h"

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// the nearest texel to the coordinates is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the nearest corresponding texel
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureNearest(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    glm::vec2 texture = texCoord;

    if (texture.x > 1) {
        texture.x = 1.0f;
    }
    if (texture.x < 0) {
        texture.x = 0.0f;
    }
    if (texture.y > 1) {
        texture.y = 1.0f;
    }
    if (texture.y < 0) {
        texture.y = 0.0f;
    }

    glm::vec2 pixelCoord = glm::vec2(floor(texture.y * image.width), floor((1.0f - texture.x) * image.height));
    int pixelIndex = image.width * (int)pixelCoord.y + (int)pixelCoord.x;

    return image.get_pixel(pixelIndex);
}

// TODO: Standard feature
// Given an image, and relevant texture coordinates, sample the texture s.t.
// a bilinearly interpolated texel is acquired from the image.
// - image;    the image object to sample from.
// - texCoord; sample coordinates, generally in [0, 1]
// - return;   the filter of the corresponding texels
// This method is unit-tested, so do not change the function signature.
glm::vec3 sampleTextureBilinear(const Image& image, const glm::vec2& texCoord)
{
    // TODO: implement this function.
    // Note: the pixels are stored in a 1D array, row-major order. You can convert from (i, j) to
    //       an index using the method seen in the lecture.
    // Note: the center of the first pixel should be at coordinates (0.5, 0.5)
    // Given texcoords, return the corresponding pixel of the image
    // The pixel are stored in a 1D array of row major order
    // you can convert from position (i,j) to an index using the method seen in the lecture
    // Note, the center of the first pixel is at image coordinates (0.5, 0.5)

    glm::vec2 texture = texCoord;

    if (texture.x > 1) {
        texture.x = 1.0f;
    }
    if (texture.x < 0) {
        texture.x = 0.0f;
    }
    if (texture.y > 1) {
        texture.y = 1.0f;
    }
    if (texture.y < 0) {
        texture.y = 0.0f;
    }

    glm::vec2 scaledTextureCoord = glm::vec2(std::max(texture.y * image.width - 0.5f, 0.0f), std::max((1.0f - texture.x) * image.height - 0.5f, 0.0f));

    scaledTextureCoord.x = std::min(scaledTextureCoord.x, image.width - 1.0f);
    scaledTextureCoord.y = std::min(scaledTextureCoord.y, image.height - 1.0f);

    glm::vec2 pixelCoordLowerLeft = glm::vec2(floor(scaledTextureCoord.x), floor(scaledTextureCoord.y));
    glm::vec2 pixelCoordUpperLeft = glm::vec2(floor(scaledTextureCoord.x), ceil(scaledTextureCoord.y));
    glm::vec2 pixelCoordLowerRight = glm::vec2(ceil(scaledTextureCoord.x), floor(scaledTextureCoord.y));
    glm::vec2 pixelCoordUpperRight = glm::vec2(ceil(scaledTextureCoord.x), ceil(scaledTextureCoord.y));

    glm::vec3 pixelLowerLeft = image.get_pixel(image.width * pixelCoordLowerLeft.y + pixelCoordLowerLeft.x);
    glm::vec3 pixelUpperLeft = image.get_pixel(image.width * pixelCoordUpperLeft.y + pixelCoordUpperLeft.x);
    glm::vec3 pixelLowerRight = image.get_pixel(image.width * pixelCoordLowerRight.y + pixelCoordLowerRight.x);
    glm::vec3 pixelUpperRight = image.get_pixel(image.width * pixelCoordUpperRight.y + pixelCoordUpperRight.x);

    float x1 = scaledTextureCoord.x - floor(scaledTextureCoord.x);
    float x2 = 1.0f - x1;

    float y1 = scaledTextureCoord.y - floor(scaledTextureCoord.y);
    float y2 = 1.0f - y1;

    glm::vec3 result = glm::vec3(pixelLowerLeft.x * x2 * y2 + pixelUpperLeft.x * x2 * y1 + pixelLowerRight.x * x1 * y2 + pixelUpperRight.x * x1 * y1,
    pixelLowerLeft.y * x2 * y2 + pixelUpperLeft.y * x2 * y1 + pixelLowerRight.y * x1 * y2 + pixelUpperRight.y * x1 * y1,
    pixelLowerLeft.z * x2 * y2 + pixelUpperLeft.z * x2 * y1 + pixelLowerRight.z * x1 * y2 + pixelUpperRight.z * x1 * y1);

    return result;
}

GLuint createTextureFromImage(std::shared_ptr<Image> image) {
    GLuint textureID;
        glGenTextures(1, &textureID);
        glBindTexture(GL_TEXTURE_2D, textureID);

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, image->width, image->height, 0,
                     (image->channels == 3 ? GL_RGB : GL_RGBA), GL_UNSIGNED_BYTE, image->get_data());

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

        glBindTexture(GL_TEXTURE_2D, 0);
        return textureID;
}

void textureVisualDebug(std::shared_ptr<Image> image, glm::vec2 hitPos, RenderState state) {
    GLuint textureID = createTextureFromImage(image);

    float width = image->width * 3;
    float height = image->height * 3;

    ImGui::Begin("Texture Mapping Visual Debug");
    ImVec2 imagePos = ImGui::GetCursorScreenPos();
    ImGui::Image(reinterpret_cast<void*>(textureID), ImVec2(width, height));

    ImDrawList* drawList = ImGui::GetWindowDrawList();

    int vertexSkip = 1;
    for (int i = 0; i < state.scene.meshes.size(); i++) {
        if (state.scene.meshes[i].vertices.size() > 1000) {
            vertexSkip = 5;
        }
        //drawList->PrimReserve(state.scene.meshes[i].triangles.size(), state.scene.meshes[i].triangles.size());
        for (int j = 0; j < state.scene.meshes[i].triangles.size(); j += vertexSkip) {

            glm::uvec3 triangle = state.scene.meshes[i].triangles[j];
            ImVec2 p1 = ImVec2(imagePos.x + state.scene.meshes[i].vertices[triangle.x].texCoord.x * width,
             imagePos.y + (1.0f - state.scene.meshes[i].vertices[triangle.x].texCoord.y) * height);
            ImVec2 p2 = ImVec2(imagePos.x + state.scene.meshes[i].vertices[triangle.y].texCoord.x * width,
             imagePos.y + (1.0f - state.scene.meshes[i].vertices[triangle.y].texCoord.y) * height);
            ImVec2 p3 = ImVec2(imagePos.x + state.scene.meshes[i].vertices[triangle.z].texCoord.x * width,
             imagePos.y + (1.0f - state.scene.meshes[i].vertices[triangle.z].texCoord.y) * height);

            drawList->AddTriangle(p1, p2, p3, IM_COL32(0, 0, 255, 255), 1.0f);
        }
    }

    drawList->AddCircleFilled(ImVec2(imagePos.x + hitPos.x * width, imagePos.y + hitPos.y * height), 5.0f, IM_COL32(255, 0, 0, 255));
    drawList->AddDrawCmd();
    ImGui::End();
    
}