#pragma once

#include "glm\glm.hpp"

#include <iostream>
#include <vector>

#include "Ray.h"

class Triangle
{
public:
	glm::vec3 color;
	glm::vec3 normal;
	std::vector<glm::vec3> vertices;

	// Constructor
	Triangle(const glm::vec3& color, const glm::vec3& normal, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& v3);

	bool rayIntersection(const Ray& ray, double& t, double& u, double& v);

	friend std::ostream& operator<<(std::ostream& out, const Triangle& triangle)
	{
		out << "V0: (" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")" << std::endl;
		out << "V1: (" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")" << std::endl;
		out << "V2: (" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")";
		return out;
	}
};