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

	Material(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_);

	glm::dvec3 brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray, const glm::dvec3& outgoing_ray);

	SurfaceType type;
	glm::dvec3 color;
	double reflection_coefficient;
};

struct Lambertian : public Material
{
	Lambertian(const SurfaceType& type_, const glm::dvec3 color_, const double& reflection_coeff_)
		: Material(type_, color_, reflection_coeff_)
	{}
};
