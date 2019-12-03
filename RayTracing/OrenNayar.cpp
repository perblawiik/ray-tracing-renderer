#include "OrenNayar.h"

OrenNayar::OrenNayar(const SurfaceType& type, const glm::dvec3 color, const double& reflection_coeff)
	: Material(type, color, reflection_coeff)
{}

OrenNayar::OrenNayar(const OrenNayar& l)
	: OrenNayar(l.surface_type, l.color, l.reflection_coefficient)
{}

glm::dvec3 OrenNayar::brdf(const glm::dvec3& surface_normal, const glm::dvec3& incoming_ray_direction, const glm::dvec3& outgoing_ray_direction)
{
	// Implementation based on Yasuhiro Fujii improvement of original Oren-Nayar formula
	// Source: https://mimosa-pudica.net/improved-oren-nayar.html
	double NdotL = glm::dot(surface_normal, outgoing_ray_direction);
	double NdotV = glm::dot(surface_normal, incoming_ray_direction);

	double s = glm::dot(outgoing_ray_direction, incoming_ray_direction) - (NdotL * NdotV);
	double t = (s > 0) ? glm::max(NdotL, NdotV) : 1.0;

	double roughness = 1.0 - reflection_coefficient;
	double denominator = glm::pi<double>() + (glm::half_pi<double>() - glm::two_thirds<double>()) * roughness;

	return color * NdotL * ((1.0 / denominator) + (roughness / denominator) * (s / t));
}