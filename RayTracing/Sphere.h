#pragma once

#include "Material.h"
#include "Ray.h"

#include "glm/glm.hpp"

#include <cmath>
#include <algorithm> 

class Sphere
{
public:
	Material* material;
	glm::dvec3 center;

	Sphere(Material* materia, const glm::dvec3& center, const double& radius);
	bool rayIntersection(const Ray& ray, double& d_near, double& d_far);

private:
	double _radius2;
};