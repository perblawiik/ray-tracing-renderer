#pragma once

#include "glm/vec3.hpp"

struct Material
{
	enum SurfaceType {
		Diffuse,
		Specular,
		Transparent,
		LightSource
	};

	// Constructor
	Material(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_);

	// Copy constructor
	Material(const Material& m);

	glm::dvec3 brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray);

	SurfaceType type;
	glm::dvec3 color;
	double reflection_coefficient;
};

struct Lambertian : public Material
{
	// Constructor
	Lambertian(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_);

	// Copy constructor
	Lambertian(const Lambertian& l);

	// Bidirectional Distribution Function
	glm::dvec3 brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray);
};
