#include "Material.h"

Material::Material(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_)
	: type(type_), color(color_), reflection_coefficient(reflection_coeff_)
{ }

glm::dvec3 Material::brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray)
{

}