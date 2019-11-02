#include "Material.h"

// Base Class Material
Material::Material(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff)
	: surface_type(type), color(color), reflection_coefficient(reflection_coeff)
{ }

Material::Material(const Material& m)
	: Material(m.surface_type, m.color, m.reflection_coefficient)
{ }


// Derived Class Lambertian Material
Lambertian::Lambertian(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff)
	: Material(type, color, reflection_coeff)
{}

Lambertian::Lambertian(const Lambertian& l)
	: Lambertian(l.surface_type, l.color, l.reflection_coefficient)
{}

glm::dvec3 Lambertian::brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray)
{
	return color / 3.14159265358979323846;
}