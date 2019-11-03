#include "Lambertian.h"

Lambertian::Lambertian(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff)
	: Material(type, color, reflection_coeff)
{}

Lambertian::Lambertian(const Lambertian& l)
	: Lambertian(l.surface_type, l.color, l.reflection_coefficient)
{}

glm::dvec3 Lambertian::brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray_direction, const glm::dvec3& outgoing_ray_direction)
{
	return color / 3.14159265358979323846;
}