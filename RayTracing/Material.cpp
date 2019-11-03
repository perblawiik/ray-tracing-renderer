#include "Material.h"

Material::Material(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff)
	: surface_type(type), color(color), reflection_coefficient(reflection_coeff)
{ }

Material::Material(const Material& m)
	: Material(m.surface_type, m.color, m.reflection_coefficient)
{ }