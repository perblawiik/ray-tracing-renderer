#pragma once

#include "Object.h"
#include "Ray.h"

#include "glm/glm.hpp"

#include <cmath>
#include <algorithm> 

class Sphere
{
public:
	glm::vec3 color;
	glm::vec3 center;

	Sphere(const glm::vec3& color, const glm::vec3& center, const double& radius);
	bool rayIntersection(const Ray& ray, double& dNear, double& dFar);

private:
	double _radius2;
};