#pragma once

#include "glm\vec3.hpp"

struct Ray
{
	glm::dvec3 start_point;
	glm::dvec3 direction;

	inline Ray(const glm::dvec3& start, const glm::dvec3& direction) : start_point(start), direction(direction) { }
};