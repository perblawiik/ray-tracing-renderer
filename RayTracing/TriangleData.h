#pragma once

#include "glm\vec3.hpp"

#include <ostream>
#include <vector>

#include "Triangle.h"

class TriangleData
{
public:
	std::vector<Triangle> triangles;

	// Constructor
	TriangleData();
	TriangleData(const double data[], const size_t num_triangles);
	void loadData(const double data[], const size_t num_triangles);
};