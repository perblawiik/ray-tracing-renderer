#include "Material.h"

// Base Class Material
Material::Material(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_)
	: type(type_), color(color_), reflection_coefficient(reflection_coeff_)
{ }

Material::Material(const Material& m)
	: Material(m.type, m.color, m.reflection_coefficient)
{ }


// Derived Class Lambertian Material
Lambertian::Lambertian(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_)
	: Material(type_, color_, reflection_coeff_)
{}

Lambertian::Lambertian(const Lambertian& l)
	: Lambertian(l.type, l.color, l.reflection_coefficient)
{}

glm::dvec3 Lambertian::brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray)
{
	return color / 3.14159265358979323846;
}