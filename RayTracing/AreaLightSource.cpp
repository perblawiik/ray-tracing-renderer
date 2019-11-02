#include "AreaLightSource.h"

AreaLightSource::AreaLightSource(const glm::dvec3& v0, const glm::dvec3& v1, const glm::dvec3& v2, Material* material, const double& watts)
	: material(material),
	intensity(watts * glm::length(glm::cross(v1 - v0, v2 - v0)) * 0.5),
	triangle(glm::normalize(glm::cross(v1 - v0, v2 - v0)), v0, v1, v2),
	_distribution(std::uniform_real_distribution<double>(0.0, 1.0))
{ }

glm::dvec3 AreaLightSource::generateRandomSamplePoint()
{
	// Barycentric coordiantes
	double u, v;

	do { // Generate random barycentric coordinates
		u = _distribution(_generator);
		v = _distribution(_generator);
	} while (u + v >= 1.0);

	// Return barycentric coordinates converted to world coordinates
	return barycentricToWorldCoordinates(u, v);
}


glm::dvec3 AreaLightSource::barycentricToWorldCoordinates(const double& u, const double& v)
{
	glm::dvec3 u_vec(u);
	glm::dvec3 v_vec(v);
	glm::dvec3 one_vec(1.0);
	return ((one_vec - u_vec - v_vec) * triangle.vertices[0]) + (u_vec * triangle.vertices[1]) + (v_vec * triangle.vertices[2]);
}