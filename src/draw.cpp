#include "draw.h"
#include <framework/opengl_includes.h>
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#ifdef __APPLE__
#include <OpenGL/GLU.h>
#else
#ifdef WIN32
// Windows.h includes a ton of stuff we don't need, this macro tells it to include less junk.
#define WIN32_LEAN_AND_MEAN
// Disable legacy macro of min/max which breaks completely valid C++ code (std::min/std::max won't work).
#define NOMINMAX
// GLU requires Windows.h on Windows :-(.
#include <Windows.h>
#endif
#include <GL/glu.h>
#endif
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/mat4x4.hpp>
DISABLE_WARNINGS_POP()
#include <algorithm>
#include "splines.h"

static void setMaterial(const Material& material)
{
    // Set the material color of the shape.
    const glm::vec4 kd4 { material.kd, 1.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, glm::value_ptr(kd4));

    const glm::vec4 zero { 0.0f };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, glm::value_ptr(zero));
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, glm::value_ptr(zero));
}

void drawExampleOfCustomVisualDebug()
{
    glBegin(GL_TRIANGLES);
    glVertex3f(0.0f, 0.0f, 0.0f);
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 1.0f, 0.0f);
    glEnd();
}


void drawTriangle (const Vertex& v0, const Vertex& v1, const Vertex& v2 ) {
    glBegin(GL_TRIANGLES);
        glNormal3fv(glm::value_ptr(v0.normal));
        glVertex3fv(glm::value_ptr(v0.position));
        glNormal3fv(glm::value_ptr(v1.normal));
        glVertex3fv(glm::value_ptr(v1.position));
        glNormal3fv(glm::value_ptr(v2.normal));
        glVertex3fv(glm::value_ptr(v2.position));
    glEnd();
}

void drawTriangle (const Vertex& v0, const Vertex& v1, const Vertex& v2, const glm::vec3& colour ) {
    glColor3f(colour.x, colour.y, colour.z);
    glBegin(GL_TRIANGLES);
        glNormal3fv(glm::value_ptr(v0.normal));
        glVertex3fv(glm::value_ptr(v0.position));
        glNormal3fv(glm::value_ptr(v1.normal));
        glVertex3fv(glm::value_ptr(v1.position));
        glNormal3fv(glm::value_ptr(v2.normal));
        glVertex3fv(glm::value_ptr(v2.position));
    glEnd();
}

void drawMesh(const Mesh& mesh)
{
    setMaterial(mesh.material);

    GLuint texture;
    
    if (mesh.material.kdTexture != nullptr) {
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        std::shared_ptr<Image> kdTexture = mesh.material.kdTexture;

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, kdTexture->width, kdTexture->height,
         0, GL_RGB, GL_UNSIGNED_BYTE, kdTexture->get_data());
        glEnable(GL_TEXTURE_2D);
    }

    glBegin(GL_TRIANGLES);
    for (const auto& triangleIndex : mesh.triangles) {
        for (int i = 0; i < 3; i++) {
            const auto& vertex = mesh.vertices[triangleIndex[i]];
            glNormal3fv(glm::value_ptr(vertex.normal));
            glTexCoord2f(vertex.texCoord.x, 1.0f - vertex.texCoord.y);
            glVertex3fv(glm::value_ptr(vertex.position));
        }
    }
    glEnd();
}

static void drawSphereInternal(const glm::vec3& center, float radius)
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    const glm::mat4 transform = glm::translate(glm::identity<glm::mat4>(), center);
    glMultMatrixf(glm::value_ptr(transform));
    auto quadric = gluNewQuadric();
    gluSphere(quadric, radius, 40, 20);
    gluDeleteQuadric(quadric);
    glPopMatrix();
}

void drawSphere(const Sphere& sphere)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    setMaterial(sphere.material);
    drawSphereInternal(sphere.center, sphere.radius);
    glPopAttrib();
}

void drawSphere(const glm::vec3& center, float radius, const glm::vec3& color /*= glm::vec3(1.0f)*/)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glColor4f(color.r, color.g, color.b, 1.0f);
    glPolygonMode(GL_FRONT, GL_FILL);
    glPolygonMode(GL_BACK, GL_FILL);
    drawSphereInternal(center, radius);
    glPopAttrib();
}

static void drawAABBInternal(const AxisAlignedBox& box)
{
    glPushMatrix();

    // front      back
    // 3 ----- 2  7 ----- 6
    // |       |  |       |
    // |       |  |       |
    // 0 ------1  4 ------5

    glBegin(GL_QUADS);
    glNormal3f(0, 0, -1);
    glVertex3f(box.lower.x, box.upper.y, box.lower.z); //3
    glVertex3f(box.upper.x, box.upper.y, box.lower.z); //2
    glVertex3f(box.upper.x, box.lower.y, box.lower.z); //1
    glVertex3f(box.lower.x, box.lower.y, box.lower.z); //0

    glNormal3f(0, 0, 1);
    glVertex3f(box.upper.x, box.lower.y, box.upper.z); //5
    glVertex3f(box.upper.x, box.upper.y, box.upper.z); //6
    glVertex3f(box.lower.x, box.upper.y, box.upper.z); //7
    glVertex3f(box.lower.x, box.lower.y, box.upper.z); //4

    glNormal3f(1, 0, 0);
    glVertex3f(box.upper.x, box.upper.y, box.lower.z); //2
    glVertex3f(box.upper.x, box.upper.y, box.upper.z); //6
    glVertex3f(box.upper.x, box.lower.y, box.upper.z); //5
    glVertex3f(box.upper.x, box.lower.y, box.lower.z); //1

    glNormal3f(-1, 0, 0);
    glVertex3f(box.lower.x, box.lower.y, box.upper.z); //4
    glVertex3f(box.lower.x, box.upper.y, box.upper.z); //7
    glVertex3f(box.lower.x, box.upper.y, box.lower.z); //3
    glVertex3f(box.lower.x, box.lower.y, box.lower.z); //0

    glNormal3f(0, 1, 0);
    glVertex3f(box.lower.x, box.upper.y, box.upper.z); //7
    glVertex3f(box.upper.x, box.upper.y, box.upper.z); //6
    glVertex3f(box.upper.x, box.upper.y, box.lower.z); //2
    glVertex3f(box.lower.x, box.upper.y, box.lower.z); //3

    glNormal3f(0, -1, 0);
    glVertex3f(box.upper.x, box.lower.y, box.lower.z); //1
    glVertex3f(box.upper.x, box.lower.y, box.upper.z); //5
    glVertex3f(box.lower.x, box.lower.y, box.upper.z); //4
    glVertex3f(box.lower.x, box.lower.y, box.lower.z); //0
    glEnd();

    glPopMatrix();
}

void drawAABB(const AxisAlignedBox& box, DrawMode drawMode, const glm::vec3& color, float transparency)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glColor4f(color.r, color.g, color.b, transparency);
    if (drawMode == DrawMode::Filled) {
        glPolygonMode(GL_FRONT, GL_FILL);
        glPolygonMode(GL_BACK, GL_FILL);
    } else {
        glPolygonMode(GL_FRONT, GL_LINE);
        glPolygonMode(GL_BACK, GL_LINE);
    }
    drawAABBInternal(box);
    glPopAttrib();
}

void drawScene(const Scene& scene)
{
    for (const auto& mesh : scene.meshes)
        drawMesh(mesh);
    for (const auto& sphere : scene.spheres)
        drawSphere(sphere);
}

void drawPlane(const glm::vec3& position, const glm::vec3& normal, const glm::vec3& up, const glm::vec3& color)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glColor4f(color.r, color.g, color.b, 0.3f); // Reduced opacity for better visibility

    // Calculate basis vectors for the plane
    glm::vec3 right = glm::normalize(glm::cross(normal, up));
    glm::vec3 planeUp = glm::normalize(glm::cross(right, normal));

    // Plane size (smaller size)
    float planeSize = 0.5f;

    glm::vec3 corners[4] = {
        position - right * planeSize - planeUp * planeSize,
        position + right * planeSize - planeUp * planeSize,
        position + right * planeSize + planeUp * planeSize,
        position - right * planeSize + planeUp * planeSize,
    };

    // Draw the plane as a quad
    glBegin(GL_QUADS);
    for (const auto& corner : corners) {
        glVertex3fv(glm::value_ptr(corner));
    }
    glEnd();

    glPopAttrib();
}

void drawEllipse(const glm::vec3& position, const glm::vec3& normal, const glm::vec3& up, float radiusX, float radiusY, const glm::vec3& color)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glColor4f(color.r, color.g, color.b, 0.5f); // Semi-transparent ellipse

    // Calculate basis vectors for the ellipse plane
    glm::vec3 right = glm::normalize(glm::cross(normal, up));
    glm::vec3 ellipseUp = glm::normalize(glm::cross(right, normal));

    glBegin(GL_LINE_LOOP); // Use GL_LINE_LOOP to draw the ellipse outline

    // Generate ellipse points
    int segments = 100; // Number of segments to approximate the ellipse
    for (int i = 0; i < segments; i++) {
        float theta = 2.0f * glm::pi<float>() * static_cast<float>(i) / static_cast<float>(segments);
        float x = radiusX * std::cos(theta);
        float y = radiusY * std::sin(theta);

        // Transform the point to world space
        glm::vec3 point = position + x * right + y * ellipseUp;
        glVertex3fv(glm::value_ptr(point));
    }

    glEnd();
    glPopAttrib();
}


void drawRay(const Ray& ray, const glm::vec3& color)
{
    const glm::vec3 hitPoint = ray.origin + std::clamp(ray.t, 0.0f, 100.0f) * ray.direction;
    const bool hit = (ray.t < std::numeric_limits<float>::max());

    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);

    glColor3fv(glm::value_ptr(color));
    glVertex3fv(glm::value_ptr(ray.origin));
    glColor3fv(glm::value_ptr(color));
    glVertex3fv(glm::value_ptr(hitPoint));
    glEnd();

    if (hit)
        drawSphere(hitPoint, 0.005f, color);

    glPopAttrib();

}

void drawLine(const glm::vec3& pos, const glm::vec3& pos2)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);

    glBegin(GL_LINES);
    glColor3f(0.5f, 0.5f, 0.5f);
    glVertex3fv(glm::value_ptr(pos));
    glVertex3fv(glm::value_ptr(pos2));
    glEnd();

    glPopAttrib();
}

void drayRayYawPitch(const glm::vec3& origin, const glm::vec2& rotation, float length) {
    float yaw = rotation.x;
    float pitch = rotation.y;

    glm::vec3 direction = glm::normalize(glm::vec3(
        cos(yaw) * cos(pitch),
        sin(pitch),
        sin(yaw) * cos(pitch)
    ));

    Ray ray;
    ray.origin = origin;
    ray.direction = direction;
    ray.t = length;

    drawRay(ray, glm::vec3(1.0f, 1.0f, 1.0f));
}

void drawSpline(const BSpline& spline, float step) {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glBegin(GL_LINE_STRIP);

    glColor3f(1, 1, 1);
    for (float t = 0; t <= 1.0f; t += step) {
        glVertex3fv(glm::value_ptr(interpolateBSpline(spline, t).position));
    }
    
    glEnd();

    for (float t = 0; t <= 1.0f; t += 5 * step) {
        SplinePoint p = interpolateBSpline(spline, t);
        drayRayYawPitch(p.position, p.rotation, 0.05f);
    }
    glPopAttrib();
}

void drawLens(const glm::vec3& position, const glm::vec3& normal, float radius, const glm::vec3& color)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glColor4f(color.r, color.g, color.b, 0.5f); // Semi-transparent lens geometry

    // Calculate basis vectors for the lens plane
    glm::vec3 right = glm::normalize(glm::cross(normal, glm::vec3(0.0f, 1.0f, 0.0f)));
    glm::vec3 lensUp = glm::normalize(glm::cross(right, normal));

    // Number of segments for drawing the lens circle
    int segments = 100;

    // Draw the lens as a filled ellipse (disk)
    glBegin(GL_TRIANGLE_FAN);
    glVertex3fv(glm::value_ptr(position)); // Center of the lens
    for (int i = 0; i <= segments; i++) {
        float theta = 2.0f * glm::pi<float>() * static_cast<float>(i) / static_cast<float>(segments);
        float x = radius * std::cos(theta);
        float y = radius * std::sin(theta);

        // Transform the point to the lens plane
        glm::vec3 point = position + x * right + y * lensUp;
        glVertex3fv(glm::value_ptr(point));
    }
    glEnd();

    glPopAttrib();
}

void drawThinLens(const glm::vec3& position, const glm::vec3& normal, float radiusX, float radiusY, const glm::vec3& color)
{
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glDisable(GL_LIGHTING);
    glColor4f(color.r, color.g, color.b, 0.5f); // Semi-transparent lens

    // Calculate basis vectors for the lens plane
    glm::vec3 right = glm::normalize(glm::cross(normal, glm::vec3(0.0f, 1.0f, 0.0f)));
    glm::vec3 lensUp = glm::normalize(glm::cross(right, normal));

    glBegin(GL_TRIANGLE_FAN);
    glVertex3fv(glm::value_ptr(position)); // Center of the lens
    int segments = 100; // Number of segments to approximate the ellipse
    for (int i = 0; i <= segments; i++) {
        float theta = 2.0f * glm::pi<float>() * static_cast<float>(i) / static_cast<float>(segments);
        float x = radiusX * std::cos(theta);
        float y = radiusY * std::sin(theta);

        // Transform the point to the lens plane
        glm::vec3 point = position + x * right + y * lensUp;
        glVertex3fv(glm::value_ptr(point));
    }
    glEnd();

    glPopAttrib();
}
