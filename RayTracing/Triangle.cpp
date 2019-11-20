#include "Triangle.h"

Triangle::Triangle(const glm::dvec3& normal, const glm::dvec3& v1, const glm::dvec3& v2, const glm::dvec3& v3)
	: normal(normal)
{
	vertices.reserve(3);
	vertices.emplace_back(v1);
	vertices.emplace_back(v2);
	vertices.emplace_back(v3);
}

bool Triangle::rayIntersection(const Ray& ray, double& t, double& u, double& v)
{
	glm::dvec3 E_1 = vertices[1] - vertices[0];
	glm::dvec3 E_2 = vertices[2] - vertices[0];
	glm::dvec3 P = glm::cross(ray.direction, E_2);

	double det = glm::dot(P, E_1);

	// Ray and triangle are parallel if length of the determinant is close to 0
	if (abs(det) < 0.0000001) {
		return false;
	}

	glm::dvec3 T = ray.start_point - vertices[0];
	double det_inversed = 1.0 / det;
	u = glm::dot(P, T) * det_inversed;

	if (u < 0.0 || u > 1.0) {
		return false;
	}

	glm::dvec3 Q = cross(T, E_1);
	v = glm::dot(Q, ray.direction) * det_inversed;

	if (v < 0.0 || u + v > 1.0) {
		return false;
	}

	t = glm::dot(Q, E_2) * det_inversed;

	// No intersections in the ray's direction if t < 0
	if (t < 0.00001) {
		return false;
	}
		
	return true;
}