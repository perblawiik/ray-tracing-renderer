#pragma once

#include "glm\glm.hpp"

#include <random>
#include <iostream>

#include "Triangle.h"

class AreaLightSource 
{
public:

	AreaLightSource(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::dvec3& color, const double& intensity);

	glm::vec3 generateRandomSamplePoint();

	// Getters
	glm::vec3 getNormal();
	glm::dvec3 getColor();

	double getIntensity();
	double getArea();

	Triangle* getTriangle();

private:
	double _intensity;
	double _area;

	std::default_random_engine _generator;
	std::uniform_real_distribution<double> _distribution;

	Triangle _triangle;

	glm::vec3 barycentricToWorldCoordinates(const double& u, const double& v);
};