#include "AreaLightSource.h"

AreaLightSource::AreaLightSource(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::dvec3& color, const double& intensity)
	: _v0(v0), _v1(v1), _v2(v2), _light_color(color), _intensity(intensity), _normal(glm::normalize(glm::cross(v1 - v0, v2 - v0)))
{ }

glm::vec3 AreaLightSource::generateRandomSamplePoint()
{
	// Barycentric coordiantes
	double u, v;

	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	// Generate random barycentric coordinates
	do {
		u = distribution(generator);
		v = distribution(generator);
	} while (u + v >= 1.0);

	// Return barycentric coordinates converted to world coordinates
	return barycentricToWorldCoordinates(u, v);
}

glm::vec3 AreaLightSource::barycentricToWorldCoordinates(const double& u, const double& v)
{
	glm::vec3 u_vec(u);
	glm::vec3 v_vec(v);
	glm::vec3 one_vec(1.0);
	return ((one_vec - u_vec - v_vec) * _v0) + (u_vec * _v1) + (v_vec * _v2);
}