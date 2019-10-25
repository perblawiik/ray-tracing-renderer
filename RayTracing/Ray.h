#pragma once

#include "glm\vec3.hpp"

struct Ray
{
	glm::vec3 start_point;
	glm::vec3 direction;

	inline Ray(const glm::vec3& start, const glm::vec3& direction) : start_point(start), direction(direction) { }
};