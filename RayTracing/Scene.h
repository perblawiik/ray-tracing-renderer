#pragma once

#include "TriangleObject.h"
#include "Sphere.h"
#include "AreaLightSource.h"

struct Scene
{
	std::vector<AreaLightSource*> light_sources;
	std::vector<TriangleObject*> triangle_objects;
	std::vector<Sphere*> spheres;

	Scene() = default;

	void addLightSource(AreaLightSource* light_source) { light_sources.emplace_back(light_source); }
	void addTriangleObject(TriangleObject* object) { triangle_objects.emplace_back(object); }
	void addSphere(Sphere* sphere) { spheres.emplace_back(sphere); }
};