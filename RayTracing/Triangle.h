#pragma once

#include "glm\glm.hpp"

#include <iostream>
#include <vector>

#include "Ray.h"

class Triangle
{
public:
	glm::dvec3 normal;
	std::vector<glm::dvec3> vertices;

	// Constructor
	Triangle(const glm::dvec3& normal, const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3);

	bool rayIntersection(const Ray& ray, double& t, double& u, double& v);

	friend std::ostream& operator<<(std::ostream& out, const Triangle& triangle)
	{
		out << "V0: (" << triangle.vertices[0].x << ", " << triangle.vertices[0].y << ", " << triangle.vertices[0].z << ")" << std::endl;
		out << "V1: (" << triangle.vertices[1].x << ", " << triangle.vertices[1].y << ", " << triangle.vertices[1].z << ")" << std::endl;
		out << "V2: (" << triangle.vertices[2].x << ", " << triangle.vertices[2].y << ", " << triangle.vertices[2].z << ")";
		return out;
	}
};