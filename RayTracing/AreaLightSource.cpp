#include "AreaLightSource.h"

AreaLightSource::AreaLightSource(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::dvec3& color, const double& intensity)
	: _intensity(intensity), _area(glm::length(glm::cross(v1 - v0, v2 - v0)) * 0.5),
	_distribution(std::uniform_real_distribution<double>(0.0, 1.0)),
	_triangle(color, glm::normalize(glm::cross(v1 - v0, v2 - v0)), v0, v1, v2)
{ 
	
}

glm::vec3 AreaLightSource::generateRandomSamplePoint()
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

Triangle* AreaLightSource::getTriangle()
{
	return &_triangle;
}

glm::vec3 AreaLightSource::getNormal()
{
	return _triangle.normal;
}

glm::dvec3 AreaLightSource::getColor()
{
	return _triangle.color;
}

double AreaLightSource::getIntensity()
{
	return _intensity;
}

double AreaLightSource::getArea()
{
	return _area;
}

glm::vec3 AreaLightSource::barycentricToWorldCoordinates(const double& u, const double& v)
{
	glm::vec3 u_vec(u);
	glm::vec3 v_vec(v);
	glm::vec3 one_vec(1.0);
	return ((one_vec - u_vec - v_vec) * _triangle.vertices[0]) + (u_vec * _triangle.vertices[1]) + (v_vec * _triangle.vertices[2]);
}