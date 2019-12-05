#pragma once

#include "Material.h"
#include "glm/gtc/type_ptr.hpp"

struct OrenNayar : public Material
{
	// Constructor
	OrenNayar(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff, const double& roughness);

	// Copy constructor
	OrenNayar(const OrenNayar& l);

	// Bidirectional Distribution Function
	glm::dvec3 brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray_direction, const glm::dvec3& outgoing_ray_direction) override;

	double roughness;
};