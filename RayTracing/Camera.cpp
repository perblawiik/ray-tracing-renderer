#include "Camera.h"

Camera::Camera(const size_t width, const size_t height, const vec3& eye)
	: _eye_point(eye), _camera_film(width, height, vec3(_eye_point.x - 1.0, 0.0, 0.0)),
	_diffuse_walls(nullptr), _specular_tetrahedron(nullptr), _specular_sphere(nullptr),
	_light_position(vec3(4.0, 0.0, 4.0)), _distribution(std::uniform_real_distribution<double>(0.0, 1.0))
{
}

void Camera::loadSceneObjects(TriangleObject* walls, Sphere* sphere, TriangleObject* tetrahedron)
{
	_diffuse_walls = walls;
	_specular_sphere = sphere;
	_specular_tetrahedron = tetrahedron;
}

void Camera::setLightPosition(const vec3& position)
{
	_light_position = position;
}

void Camera::addLightSource(AreaLightSource* light_source) 
{
	_light_source = light_source;
}

void Camera::render()
{
	for (int j = 0; j < _camera_film.height; ++j) {
		// Print progress for every full percent
		if (j % 8 == 0) {
			float progress = ((float)j / (float)_camera_film.height) * 100.0;
			std::cout << "Progress: " << progress << " %" << std::endl;
		}

		for (int i = 0; i < _camera_film.width; ++i) {
			// Compute current pixel coordinates to world coordinates
			vec3 pixel(_camera_film.position.x, (double)i * 0.0025 - 0.99875, -((double)j * 0.0025 - 0.99875));

			// Define the ray from camera eye to pixel
			Ray ray(_eye_point, normalize(pixel - _eye_point));

			// Number of ray reflections to trace
			const int max_reflections = 4;

			// Number of hemispherical sample rays
			const int N = 64;

			// Compute the incoming radiance
			dvec3 final_color = tracePath(ray, max_reflections, N);

			// Save pixel color
			_camera_film.pixel_data.emplace_back(final_color);
		}
	}
}

void Camera::createImage(const char* file_name)
{
	size_t image_width = _camera_film.width;
	size_t image_height = _camera_film.height;

	// Allocate memory for the bitmap file
	Bitmap image(image_width, image_height);

	int index = 0;
	for (int j = 0; j < image_height; ++j) {
		for (int i = 0; i < image_width; ++i) {
			// Fetch current pixel intensity
			dvec3 pixel = _camera_film.pixel_data[index++];

			// Convert color values to unsigned char
			uint8_t r = uint8_t(min(pixel.x, 1.0) * 255);
			uint8_t g = uint8_t(min(pixel.y, 1.0) * 255);
			uint8_t b = uint8_t(min(pixel.z, 1.0) * 255);

			// Set pixel color in the bitmap image
			image.setPixelColor(i, j, r, g, b, 255);
		}
	}

	// Write data to a bitmap file
	image.writeToFile(file_name);
}

//***************** PRIVATE ****************//

dvec3 Camera::tracePath(const Ray& ray, const int reflection_count, const int mc_sample_ray_count)
{
	// The closest intersection for current ray
	IntersectionPoint closest_point;

	// Test the ray against all triangle objects
	triangleIntersectionTests(ray, closest_point);

	// Test the ray against the sphere
	sphereIntersectionTest(ray, closest_point);

	// If distance is negative, no intersection was found in the ray direction
	if (closest_point.distance < 0.0) {
		return dvec3(0.0);
	}
	else {
		// Light source
		if (closest_point.type == IntersectionPoint::SurfaceType::LightSource) {
			return _light_source->getColor();
		}
		// Specular surface
		else if (closest_point.type == IntersectionPoint::SurfaceType::Specular) {
			// Base case (if reflection count is zero end recursion)
			if (reflection_count <= 0)
				return dvec3(0.0);

			// Compute reflected ray
			vec3 reflect_direction = reflect(ray.direction, closest_point.normal);
			Ray reflected_ray(closest_point.position, reflect_direction);

			double cos_theta = max(dot(closest_point.normal, reflected_ray.direction), 0.0);

			// Recursive path tracing (perfect reflection)
			return tracePath(reflected_ray, reflection_count - 1, mc_sample_ray_count) * cos_theta;
		} 
		else { // Diffuse surface
			// Lambertian BRDF
			dvec3 brdf = closest_point.color / PI;

			// Compute direct light from light source
			dvec3 direct_light = computeDirectLight(closest_point, brdf, 20);

			// Base case (if reflection count is zero, only return the direct light)
			if (reflection_count <= 0)
				return direct_light;

			vec3 local_x_axis(0.0);
			vec3 local_z_axis(0.0);
			// Create a local coordinate system for the hemisphere of the intersection point with the surface normal as the y-axis
			createLocalCoordinateSystem(closest_point.normal, local_x_axis, local_z_axis);

			// Perform a Monte Carlo integration for indirect lighting with N samples
			dvec3 indirect_light(0.0);
			for (size_t n = 0; n < mc_sample_ray_count; ++n) {
				// Generate and compute a random direction in the hemisphere
				double cos_theta = _distribution(_generator); // Let the first random number be equal to cos(theta)
				double random_2 = _distribution(_generator);

				vec3 sample_direction_local = hemisphereSampleDirection(cos_theta, random_2);

				// Transform direction from local to world by multiplying with the local coordinate system matrix
				// Note that we only transform the direction so the translation part is not needed (in that case we would use a 4x4 matrix)
				vec3 sample_direction_world(
					sample_direction_local.x * local_z_axis.x + sample_direction_local.y * closest_point.normal.x + sample_direction_local.z * local_x_axis.x,
					sample_direction_local.x * local_z_axis.y + sample_direction_local.y * closest_point.normal.y + sample_direction_local.z * local_x_axis.y,
					sample_direction_local.x * local_z_axis.z + sample_direction_local.y * closest_point.normal.z + sample_direction_local.z * local_x_axis.z
				);

				// Create a ray from the sample direction
				Ray sample_ray(closest_point.position, sample_direction_world);

				// Note that the first random variable corresponds to cos(theta) which gives us the diffuse distribution
				indirect_light += tracePath(sample_ray, 0, mc_sample_ray_count / 2) * cos_theta;
			}

			// The indirect light needs to be divided by the pdf constant (probability density function)
			// Note that the pdf is constant in this case because all of the random directions share the same probability (equiprobability)
			double pdf = 1.0 / (2.0 * PI);
			// Also divide the sum by N to complete the Monte Carlo sampling
			indirect_light /= ((double)mc_sample_ray_count * pdf);

			// Multiply the brdf surface color with the light
			indirect_light *= brdf;

			// Combine the direct and indirect light and multiply with surface color
			return (direct_light + indirect_light);
		}
	}
}

bool Camera::shadowRay(const vec3& surface_point, const vec3& point_to_light, const double& light_distance)
{
	// Define a ray pointing towards the light source
	Ray light_ray(surface_point, point_to_light);
	IntersectionPoint light_path_point;

	bool intersection_found = false;

	for (Triangle triangle : _specular_tetrahedron->triangles) {
		double t, u, v = -1.0;
		intersection_found = triangle.rayIntersection(light_ray, t, u, v);
		if (intersection_found) {
			light_path_point.distance = t;
			break;
		}
	}

	// Shadow
	if (intersection_found && light_path_point.distance < light_distance) {
		return true;
	}

	double d_near, d_far = -1.0;
	intersection_found = _specular_sphere->rayIntersection(light_ray, d_near, d_far);
	if (intersection_found && d_near < light_path_point.distance) {
		light_path_point.distance = d_near;
	}

	// Shadow
	if (intersection_found && light_path_point.distance < light_distance) {
		return true;
	}

	// No shadow
	return false;
}

dvec3 Camera::computeDirectLight(const IntersectionPoint& surface_point, const dvec3& brdf, const size_t sample_ray_count)
{
	dvec3 direct_light(0.0);
	std::vector<vec3> light_sample_points;

	// Generate sample points on the area light source
	for (size_t i = 0; i < sample_ray_count; i++) {
		light_sample_points.push_back(_light_source->generateRandomSamplePoint());
	}

	// Cast a shadow ray for each sample point
	for (auto sample_point : light_sample_points) {
		vec3 point_to_light_direction = normalize(sample_point - surface_point.position);
		double point_to_light_distance = distance(sample_point, surface_point.position);

		if (!shadowRay(surface_point.position, point_to_light_direction, point_to_light_distance)) {
			double cos_theta_out = max(dot(surface_point.normal, point_to_light_direction), 0.0);
			double cos_theta_in = dot(-point_to_light_direction, _light_source->getNormal());
			
			direct_light += brdf * cos_theta_out * cos_theta_in / (point_to_light_distance * point_to_light_distance);
		}
	}

	return (direct_light * _light_source->getColor() * _light_source->getArea() * _light_source->getIntensity()) / (double) sample_ray_count;
}

void Camera::createLocalCoordinateSystem(const vec3& normal, vec3& local_x_axis, vec3& local_z_axis)
{
	// If the normals y-coordinate is smaller than the x-coordinate, next axis should lie in the y-plane (y = 0)
	if (abs(normal.x) > abs(normal.y)) {
		local_x_axis = vec3(normal.z, 0.0, -normal.x) / sqrt(normal.x * normal.x + normal.z * normal.z);
	}
	else {
		local_x_axis = vec3(0.0, -normal.z, normal.y) / sqrt(normal.y * normal.y + normal.z * normal.z);
	}

	local_z_axis = cross(normal, local_x_axis);
}

vec3 Camera::hemisphereSampleDirection(const double &cos_theta, const double &random_2)
{
	// We assume that the first random value is cos(theta), which is equal to the y-coordiante
	// Theta is the inclination angle and phi is the azimuth angle
	double sin_theta = sqrt(1 - cos_theta * cos_theta);
	double phi = 2 * PI * random_2;

	float x = sin_theta * cos(phi);
	float z = sin_theta * sin(phi);

	return vec3(x, cos_theta, z);
}

void Camera::triangleIntersectionTests(const Ray& ray, IntersectionPoint& closest_point)
{
	// Walls, ceiling and floor
	for (Triangle triangle : _diffuse_walls->triangles) {
		double t, u, v = -1.0;
		if (triangle.rayIntersection(ray, t, u, v) == true) {
			// Depth test (if we have multiple intersections, save the closest)
			if (closest_point.distance < 0.0 || closest_point.distance > t) {
				// Save closest distance
				closest_point.distance = t;
				// Compute world coordinates from barycentric triangle coordinates
				closest_point.position = barycentricToWorldCoordinates(triangle, u, v);
				closest_point.normal = triangle.normal;
				closest_point.color = triangle.color;
				closest_point.type = IntersectionPoint::SurfaceType::Diffuse;
			}
		}
	}
	
	// Tetrahedron
	for (Triangle triangle : _specular_tetrahedron->triangles) {
		double t, u, v = -1.0;
		if (triangle.rayIntersection(ray, t, u, v) == true) {
			// Depth test (if we have multiple intersections, save the closest)
			if (closest_point.distance < 0.0 || closest_point.distance > t) {
				// Save closest distance
				closest_point.distance = t;
				// Compute world coordinates from barycentric triangle coordinates
				closest_point.position = barycentricToWorldCoordinates(triangle, u, v);
				closest_point.normal = triangle.normal;
				closest_point.color = triangle.color;
				closest_point.type = IntersectionPoint::SurfaceType::Specular;
			}
		}
	}

	// Area Light Source
	double t, u, v = -1.0;
	auto lightTriangle = _light_source->getTriangle();
	if (lightTriangle->rayIntersection(ray, t, u, v) == true) {
		// Depth test (if we have multiple intersections, save the closest)
		if (closest_point.distance < 0.0 || closest_point.distance > t) {
			// Save closest distance
			closest_point.distance = t;
			// Compute world coordinates from barycentric triangle coordinates
			closest_point.position = barycentricToWorldCoordinates(*lightTriangle, u, v);
			closest_point.normal = lightTriangle->normal;
			closest_point.color = lightTriangle->color;
			closest_point.type = IntersectionPoint::SurfaceType::LightSource;
		}
	}
}

void Camera::sphereIntersectionTest(const Ray& ray, IntersectionPoint& closest_point)
{
	double d_near, d_far = -1.0;
	if (_specular_sphere->rayIntersection(ray, d_near, d_far) == true) {
		// Depth test (if we have multiple intersections, save the closest)
		if (closest_point.distance < 0.0 || closest_point.distance > d_near) {
			// Save closest distance
			closest_point.distance = d_near;
			// Compute the intersection point position
			closest_point.position = ray.start_point + vec3(d_near) * ray.direction;

			// Set attributes
			closest_point.normal = normalize(closest_point.position - _specular_sphere->center);
			closest_point.color = _specular_sphere->color;
			closest_point.type = IntersectionPoint::SurfaceType::Specular;
		}
	}
}

vec3 Camera::barycentricToWorldCoordinates(const Triangle& triangle, const double& u, const double& v) 
{
	vec3 u_vec(u);
	vec3 v_vec(v);
	vec3 one_vec(1.0);
	return ((one_vec - u_vec - v_vec) * triangle.vertices[0]) + (u_vec * triangle.vertices[1]) + (v_vec * triangle.vertices[2]);
}

double Camera::max(const double& a, const double& b) 
{
	return (a < b) ? b : a;
}