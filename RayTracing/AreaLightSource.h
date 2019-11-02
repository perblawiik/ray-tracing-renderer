#pragma once

#include "glm\glm.hpp"

#include <random>
#include <iostream>

#include "Triangle.h"
#include "Material.h"

class AreaLightSource 
{
public:

	AreaLightSource(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, Material* material, const double& watts);

	glm::dvec3 generateRandomSamplePoint();

	Material* material;
	double intensity;
	Triangle triangle;

private:

	std::default_random_engine _generator;
	std::uniform_real_distribution<double> _distribution;

	glm::dvec3 barycentricToWorldCoordinates(const double& u, const double& v);
};