#pragma once

#include "glm\glm.hpp"

#include <random>

class AreaLightSource 
{
public:

	AreaLightSource(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::dvec3& color, const double& intensity);

	glm::vec3 generateRandomSamplePoint();

private:
	glm::vec3 _v0;
	glm::vec3 _v1;
	glm::vec3 _v2;

	glm::dvec3 _light_color;
	double _intensity;

	glm::vec3 _normal;

	glm::vec3 barycentricToWorldCoordinates(const double& u, const double& v);
};