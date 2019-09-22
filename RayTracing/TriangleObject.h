#pragma once

#include "glm\vec3.hpp"

#include <ostream>
#include <vector>

#include "Triangle.h"

class TriangleObject
{
public:
	std::vector<Triangle> triangles;

	// Constructor
	TriangleObject();
	TriangleObject(const double data[], const size_t num_triangles);
	void loadData(const double data[], const size_t num_triangles);
	void createTetrahedron(const glm::vec3& color, const glm::vec3& position, const double& scale);
};