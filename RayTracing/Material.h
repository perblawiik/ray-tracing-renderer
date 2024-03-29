#pragma once

#include "glm/vec3.hpp"
#include "glm/gtc/constants.hpp"

struct Material
{
	enum SurfaceType {
		Diffuse,
		Specular,
		Transparent,
		LightSource
	};

	// Constructor
	Material(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff);

	// Copy constructor
	Material(const Material& m);

	// Destructor
	virtual ~Material() = default;

	virtual glm::dvec3 brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray) = 0;

	SurfaceType surface_type;
	glm::dvec3 color;
	double reflection_coefficient;
};