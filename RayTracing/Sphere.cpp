#include "Sphere.h"

Sphere::Sphere(const vec3& color, const vec3& center, const double& radius)
	: color(color), center(center), _radius(radius)
{}

bool Sphere::rayIntersection(const Ray& ray, double& d_near, double& d_far)
{
	vec3 center_to_ray_origin = ray.start_point - center;

	double b = dot(ray.direction, center_to_ray_origin);
	double c = dot(center_to_ray_origin, center_to_ray_origin) - (_radius * _radius);
	double sqr_bb_min_c = sqrt(b*b - c);

	d_near = -b + sqr_bb_min_c;
	d_far = -b - sqr_bb_min_c;

	if (d_near > d_far)
		std::swap(d_near, d_far);

	if (d_near < 0.0)
		return false;

	return true;
}