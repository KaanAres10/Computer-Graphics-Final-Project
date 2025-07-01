#include "splines.h"

SplinePoint interpolateBSpline(const BSpline& bSpline, float t)
{
    if (!bSpline.repeatKnots) {
        switch (bSpline.type) {
        case Linear:
            return SplinePoint { interpolateUniformLinear(bSpline.position, t), interpolateUniformLinear(bSpline.rotation, t) };
        case Quadratic:
            return SplinePoint { interpolateUniformQuadratic(bSpline.position, t), interpolateUniformQuadratic(bSpline.rotation, t) };
        case Cubic:
            return SplinePoint { interpolateUniformCubic(bSpline.position, t), interpolateUniformCubic(bSpline.rotation, t) };
        default:
            return { { 0, 0, 0 }, { 0, 0 } };
        }
    } else {
        switch (bSpline.type) {
        case Linear:
            return SplinePoint { interpolateUniformLinear(bSpline.position, t), interpolateUniformLinear(bSpline.rotation, t) };
        case Quadratic:
            return SplinePoint { interpolateConnectedQuadratic(bSpline.position, t), interpolateConnectedQuadratic(bSpline.rotation, t) };
        case Cubic:
            return SplinePoint { interpolateConnectedCubic(bSpline.position, t), interpolateConnectedCubic(bSpline.rotation, t) };
        default:
            return { { 0, 0, 0 }, { 0, 0 } };
        }
    }
}

template <typename T>
T interpolateUniformLinear(const std::vector<T>& controlPoints, float t)
{
    t *= float(controlPoints.size() - 1);
    t += 1;
    T sum { 0 };

    for (int i = 0; i < controlPoints.size(); i++) {
        T pos = controlPoints[i];
        if (t >= i && t < i + 1) {
            sum += pos * (t - i);
        }

        else if (i + 1 <= t && t <= i + 2) {
            sum += pos * (2 - t + i);
        }
    }

    return sum;
}

template <typename T>
T interpolateUniformQuadratic(const std::vector<T>& controlPoints, float t)
{
    t *= float(controlPoints.size() - 2);
    t += 2;
    T sum { 0 };

    for (int i = 0; i < controlPoints.size(); i++) {
        T pos = controlPoints[i];
        if (t >= i && t < i + 1) {
            sum += (pos * ((t - i)*(t-i) * 0.5f));
        }

        else if (i + 1 <= t && t <= i + 2) {
            float u = t - i - 1;
            u = -u * u + u + 0.5f;
            sum += pos * u;
        } 
        
        else if (i+2<=t&&t<i+3) {
            float u = t - i - 2;
            u = (1 - u) * (1 - u);
            u *= 0.5f;
            sum += pos * u;
        }
    }

    return sum;
}

template <typename T>
T interpolateConnectedQuadratic(const std::vector<T>& controlPoints, float t)
{
    T sum(0);
    t *= float(controlPoints.size() + 2);
    
    if (t < 1) {
        sum += controlPoints[0] * (1.0f - t) * (1.0f - t);
        sum += controlPoints[1] * (2.0f * t - 3.0f / 2.0f * t * t);
    } else if (t < 2) {
        sum += controlPoints[1] * (1.0f - t + 1) * (1.0f - t + 1) / 2.0f;
    }

    if (t > float(controlPoints.size())) {
        float t2 = t - float(controlPoints.size());
        t2 = 2.0f - t2;

        if (t2 < 1) {
            sum += controlPoints[controlPoints.size()-1] * (1.0f - t2) * (1.0f - t2);
            sum += controlPoints[controlPoints.size()-2] * (2.0f * t2 - 3.0f / 2.0f * t2 * t2);
        } else if (t2 < 2) {
            sum += controlPoints[controlPoints.size()-2] * (1.0f - t2 + 1) * (1.0f - t2 + 1) / 2.0f;
        }
    }

    for (int i = 0; i < controlPoints.size(); i++) {
        T pos = controlPoints[i];
        if (t >= i && t < i + 1) {
            sum += (pos * ((t - i) * (t - i) * 0.5f));
        }

        else if (i + 1 <= t && t <= i + 2) {
            float u = t - i - 1;
            u = -u * u + u + 0.5f;
            sum += pos * u;
        }

        else if (i + 2 <= t && t < i + 3) {
            float u = t - i - 2;
            u = (1 - u) * (1 - u);
            u *= 0.5f;
            sum += pos * u;
        }
    }

    return sum;
}

template <typename T>
T interpolateUniformCubic(const std::vector<T>& controlPoints, float t)
{
    t *= float(controlPoints.size() - 3);
    t += 3;
    T sum { 0 };

    for (int i = 0; i < controlPoints.size(); i++) {
        T pos = controlPoints[i];
        if (t >= i && t < i + 1) {
            sum += (pos * ((t - i) * (t - i) * (t - i) / 6.0f));
        }

        else if (i + 1 <= t && t <= i + 2) {
            float u = t - i - 1;
            u = -3.0 * u * u * u + 3 * u * u + 3 * u + 1.0f;
            sum += pos * u / 6.0f;
        }

        else if (i + 2 <= t && t < i + 3) {
            float u = t - i - 2;
            u = 3 * u * u * u - 6 * u * u + 4;
            sum += pos * u / 6.0f;
        }

        else if (i + 3 <= t && t < i + 4) {
            float u = t - i - 3;
            u = -u * u * u + 3 * u * u - 3 * u + 1;

            sum += pos * u / 6.0f;
        }
    }

    return sum;
}

template <typename T>
T interpolateConnectedCubic(const std::vector<T>& controlPoints, float t)
{   
    return interpolateUniformCubic(controlPoints, t);
}

