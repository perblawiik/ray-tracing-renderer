#pragma once

#include "glm\vec3.hpp"

#include <ostream>
#include <vector>

#include "Triangle.h"
#include "Material.h"

class TriangleObject
{
public:
	std::vector<Triangle> triangles;
	Material* material;

	// Constructor
	TriangleObject(Material *material);

	void loadData(const double data[], const size_t num_triangles);
	void createTetrahedron(const glm::dvec3& position, const double& scale);
};