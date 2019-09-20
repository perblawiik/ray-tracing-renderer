#pragma once

#include <vector>

#include "TriangleData.h"
#include "Sphere.h"

class Scene
{
public:
	std::vector<TriangleData> triangles;
	std::vector<Sphere> spheres;
private:
};