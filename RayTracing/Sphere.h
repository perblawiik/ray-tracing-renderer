#pragma once

#include "Ray.h"

#include "glm/glm.hpp"

#include <cmath>
#include <algorithm> 
#include <iostream>

using namespace glm;

class Sphere
{
public:
	vec3 color;
	vec3 center;

	Sphere(const vec3& color, const vec3& center, const double& radius);
	bool rayIntersection(const Ray& ray, double& dNear, double& dFar);

private:
	double _radius;
};