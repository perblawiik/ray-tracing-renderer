#pragma once

#include "Material.h"

struct Lambertian : public Material
{
	// Constructor
	Lambertian(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff);

	// Copy constructor
	Lambertian(const Lambertian& l);

	// Bidirectional Distribution Function
	glm::dvec3 brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray_direction, const glm::dvec3& outgoing_ray_direction) override;
};