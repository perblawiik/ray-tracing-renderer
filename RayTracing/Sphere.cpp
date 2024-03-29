#include "Sphere.h"

Sphere::Sphere(Material* material, const glm::dvec3& center, const double& radius)
	: material(material), center(center), _radius2(radius * radius)
{}

bool Sphere::rayIntersection(const Ray& ray, double& d_near, double& d_far)
{
	glm::dvec3 ray_origin_to_center = center - ray.start_point;
	double cdir_dot_rdir = glm::dot(ray_origin_to_center, ray.direction);

	// ray to center vector and ray direction is facing opposite directions
	if (cdir_dot_rdir < 0.0)
		return false;

	// Compute the square of center_to_ray_dist directly since we dont need the square root for anything
	double center_to_ray_dist2 = glm::dot(ray_origin_to_center, ray_origin_to_center) - (cdir_dot_rdir * cdir_dot_rdir);

	// Missed intersection (distance to ray is larger than the radius of the sphere)
	if (center_to_ray_dist2 > _radius2)
		return false;

	// The distance from intersection point to the center in the ray direction
	double hit_to_center_dist = sqrt(_radius2 - center_to_ray_dist2);

	d_near = cdir_dot_rdir - hit_to_center_dist;
	d_far = cdir_dot_rdir + hit_to_center_dist;

	if (d_near > d_far)
		std::swap(d_near, d_far);

	// If d_near is negative, ray origin might be inside the sphere
	if (d_near < 0.0) {
		d_near = d_far;

		// If closest positive intersection point is negative, the intersection is before ray origin
		if (d_near < 0.0)
			return false;
	}

	return true;
}