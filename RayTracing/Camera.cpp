#include "Camera.h"

Camera::Camera(const size_t width, const size_t height, const dvec3& eye)
	: _eye_point(eye), _camera_film(width, height, dvec3(_eye_point.x - 1.0, 0.0, 0.0)), _scene(nullptr),
	_distribution(std::uniform_real_distribution<double>(0.0, 1.0))
{
	setupCameraMatrix();
	computePixelWidth();
}

void Camera::loadScene(Scene* scene)
{
	_scene = scene;
}

void Camera::render(const int& num_samples)
{
	// Image dimensions
	const size_t height = _camera_film.height;
	const size_t width = _camera_film.width;

	// Number of ray reflections to trace
	const int max_reflections = 5;

	// Samples per pixel
	for (int n = 0; n < num_samples; ++n) {
		// Print out the progress in percentage
		float progress = (float)n / (float)num_samples;
		std::cout << "Progress: " << progress * 100.0 << " %" << std::endl;

		// Parallization with OpenMP
		#pragma omp parallel for schedule(dynamic)
		for (int pixel_y = 0; pixel_y < height; ++pixel_y) {
			for (int pixel_x = 0; pixel_x < width; ++pixel_x) {
				// Normalize coordinates by transforming from [0 : width/height] to the range [-1.0 : 1.0]
				dvec2 pixel_normalized = normalizedPixelCoord(pixel_x, pixel_y);

				// Generate a random offset in the pixel plane (used for ray randomization)
				double dx = _distribution(_generator) * _pixel_width - (0.5 * _pixel_width);
				double dy = _distribution(_generator) * _pixel_width - (0.5 * _pixel_width);

				// Replace aliasing with noise by adding the random offset to the pixel
				dvec4 pixel_sample(pixel_normalized.x + dx, pixel_normalized.y + dy, 1.0, 1.0);

				// Transform screen coordinates to world coordinates
				dvec3 pixel_world(_transform_matrix * pixel_sample);

				// Define the ray from camera eye to pixel
				Ray ray(_eye_point, normalize(pixel_world));

				// Compute the incoming radiance
				dvec3 final_color = tracePath(ray, max_reflections);

				// Clamp the values between [0.0, 1.0]
				clamp(final_color, 0.0, 1.0);
	
				// Save pixel color
				_camera_film.addPixelData(pixel_x, pixel_y, final_color / (double)num_samples);
			}
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

dvec3 Camera::tracePath(const Ray& ray, const int reflection_count)
{
	// Base case (if reflection count is zero end recursion)
	if (reflection_count == 0)
		return dvec3(0.0);

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

	// Light source
	if (closest_point.material->surface_type == Material::SurfaceType::LightSource) {
		return closest_point.material->color;
	}
	// Specular surface
	else if (closest_point.material->surface_type == Material::SurfaceType::Specular) {
		// Compute reflected ray
		dvec3 reflect_direction = reflect(ray.direction, closest_point.normal);
		Ray reflected_ray(closest_point.position, reflect_direction);

		// Recursive path tracing (perfect reflection)
		return tracePath(reflected_ray, reflection_count - 1);
	} 
	// Diffuse surface
	else if (closest_point.material->surface_type == Material::SurfaceType::Diffuse) {
		// Let the first random number be equal to cos(theta)
		double cos_theta = _distribution(_generator); 

		// Russian roulette
		double phi = (cos_theta * 2.0 * PI) / closest_point.material->reflection_coefficient;
		if (phi > 2.0 * PI) {
			return dvec3(0.0);
		}

		double random_2 = _distribution(_generator);

		// Generate and compute a random direction in the hemisphere
		dvec3 sample_direction_world = hemisphereSampleDirection(cos_theta, random_2, closest_point.normal);

		// Create a ray from the sample direction
		Ray sample_ray(closest_point.position, sample_direction_world);

		// BRDF
		dvec3 brdf = closest_point.material->brdf(closest_point.normal, ray.direction, sample_ray.direction);

		// Recursion
		dvec3 indirect_light = tracePath(sample_ray, reflection_count - 1) * brdf * phi;
		
		// The indirect light needs to be divided by the pdf constant (probability density function)
		// Note that the pdf is constant in this case because all of the random directions share the same probability (equiprobability)
		//double pdf = 1.0 / (2.0 * PI);
		// Also divide the sum by N to complete the Monte Carlo sampling
		//indirect_light /= pdf;

		// Compute direct light from light source
		dvec3 direct_light = computeDirectLight(closest_point, brdf, 20);

		// Combine the direct and indirect light and multiply with surface color
		return (direct_light + indirect_light);
	}
	else {
		return dvec3(0.0);
	}
}

bool Camera::shadowRay(const dvec3& surface_point, const dvec3& point_to_light, const double& light_distance)
{
	// Define a ray pointing towards the light source
	Ray light_ray(surface_point, point_to_light);

	// Only check shadow ray on the tetrahedron (last item in the triangle_objects array)
	auto tetrahedron = _scene->triangle_objects.back();
	for (Triangle triangle : tetrahedron->triangles) {
		double t, u, v = -1.0;
		bool intersection_found = triangle.rayIntersection(light_ray, t, u, v);

		// Shadow
		if (intersection_found && t < light_distance) {
			return true;
		}
	}

	for (Sphere* sphere : _scene->spheres) {
		double d_near, d_far = -1.0;
		bool intersection_found = sphere->rayIntersection(light_ray, d_near, d_far);

		// Shadow
		if (intersection_found && d_near < light_distance) {
			return true;
		}
	}

	// No shadow
	return false;
}

dvec3 Camera::computeDirectLight(const IntersectionPoint& surface_point, const dvec3& brdf, const size_t sample_ray_count)
{
	dvec3 total_light(0.0);

	// Sum the direct light from all light sources
	for (auto light_source : _scene->light_sources) {
		dvec3 direct_light(0.0);
		std::vector<dvec3> light_sample_points;

		// Generate sample points on the area light
		for (size_t i = 0; i < sample_ray_count; i++) {
			light_sample_points.push_back(light_source->generateRandomSamplePoint());
		}

		// Cast a shadow ray for each sample point
		for (auto sample_point : light_sample_points) {
			dvec3 point_to_light_direction = normalize(sample_point - surface_point.position);
			double point_to_light_distance = distance(sample_point, surface_point.position);

			if (!shadowRay(surface_point.position, point_to_light_direction, point_to_light_distance)) {
				double cos_theta_out = max(dot(surface_point.normal, point_to_light_direction), 0.0);
				double cos_theta_in = dot(-point_to_light_direction, light_source->triangle.normal);

				direct_light += brdf * cos_theta_out * cos_theta_in / (point_to_light_distance * point_to_light_distance);
			}
		}
		// Average the light
		direct_light *= (light_source->material->color * light_source->intensity / (double)sample_ray_count);
		total_light += direct_light;
	}
	
	return total_light;
}

void Camera::createLocalCoordinateSystem(const dvec3& normal, dvec3& local_x_axis, dvec3& local_z_axis)
{
	// If the normals y-coordinate is smaller than the x-coordinate, next axis should lie in the y-plane (y = 0)
	if (abs(normal.x) > abs(normal.y)) {
		local_x_axis = dvec3(normal.z, 0.0, -normal.x) / sqrt(normal.x * normal.x + normal.z * normal.z);
	}
	else {
		local_x_axis = dvec3(0.0, -normal.z, normal.y) / sqrt(normal.y * normal.y + normal.z * normal.z);
	}

	local_z_axis = cross(normal, local_x_axis);
}

dvec3 Camera::hemisphereSampleDirection(const double &cos_theta, const double &random_2, const dvec3& surface_normal)
{
	// We assume that the first random value is cos(theta), which is equal to the y-coordiante
	// Theta is the inclination angle and phi is the azimuth angle
	double sin_theta = sqrt(1 - cos_theta * cos_theta);
	double phi = 2 * PI * random_2;

	float x = sin_theta * cos(phi);
	float z = sin_theta * sin(phi);

	dvec3 sample_direction_local(x, cos_theta, z);

	// Create a local coordinate system for the hemisphere of the intersection point with the surface normal as the y-axis
	dvec3 local_x_axis(0.0);
	dvec3 local_z_axis(0.0);
	createLocalCoordinateSystem(surface_normal, local_x_axis, local_z_axis);

	// Transform direction from local to world by multiplying with the local coordinate system matrix
	// Note that we only transform the direction so the translation part is not needed (in that case we would use a 4x4 matrix)
	dvec3 sample_direction_world(
		sample_direction_local.x * local_z_axis.x + sample_direction_local.y * surface_normal.x + sample_direction_local.z * local_x_axis.x,
		sample_direction_local.x * local_z_axis.y + sample_direction_local.y * surface_normal.y + sample_direction_local.z * local_x_axis.y,
		sample_direction_local.x * local_z_axis.z + sample_direction_local.y * surface_normal.z + sample_direction_local.z * local_x_axis.z
	);

	return sample_direction_world;
}

void Camera::triangleIntersectionTests(const Ray& ray, IntersectionPoint& closest_point)
{
	// Check triangle objects
	for (TriangleObject* object : _scene->triangle_objects) {
		for (Triangle triangle : object->triangles) {
			double t, u, v = -1.0;
			if (triangle.rayIntersection(ray, t, u, v) == true) {
				// Depth test (if we have multiple intersections, save the closest)
				if (closest_point.distance < 0.0 || closest_point.distance > t) {
					// Save closest distance
					closest_point.distance = t;
					// Compute world coordinates from barycentric triangle coordinates
					closest_point.position = barycentricToWorldCoordinates(triangle, u, v);
					closest_point.normal = triangle.normal;
					closest_point.material = object->material;
				}
			}
		}
	}

	// Area Light Source
	for (auto light_source : _scene->light_sources) {
		double t, u, v = -1.0;
		if (light_source->triangle.rayIntersection(ray, t, u, v) == true) {
			// Depth test (if we have multiple intersections, save the closest)
			if (closest_point.distance < 0.0 || closest_point.distance > t) {
				// Save closest distance
				closest_point.distance = t;
				// Compute world coordinates from barycentric triangle coordinates
				closest_point.position = barycentricToWorldCoordinates(light_source->triangle, u, v);
				closest_point.normal = light_source->triangle.normal;
				closest_point.material = light_source->material;
			}
		}
	}
}

void Camera::sphereIntersectionTest(const Ray& ray, IntersectionPoint& closest_point)
{
	for (Sphere* sphere : _scene->spheres) {
		double d_near, d_far = -1.0;
		if (sphere->rayIntersection(ray, d_near, d_far) == true) {
			// Depth test (if we have multiple intersections, save the closest)
			if (closest_point.distance < 0.0 || closest_point.distance > d_near) {
				// Save closest distance
				closest_point.distance = d_near;
				// Compute the intersection point world coordinate position
				closest_point.position = ray.start_point + dvec3(d_near) * ray.direction;

				// Set attributes
				closest_point.normal = normalize(closest_point.position - sphere->center);
				closest_point.material = sphere->material;
			}
		}
	}
}

dvec3 Camera::barycentricToWorldCoordinates(const Triangle& triangle, const double& u, const double& v) 
{
	dvec3 u_vec(u);
	dvec3 v_vec(v);
	dvec3 one_vec(1.0);
	return ((one_vec - u_vec - v_vec) * triangle.vertices[0]) + (u_vec * triangle.vertices[1]) + (v_vec * triangle.vertices[2]);
}

double Camera::max(const double& a, const double& b) 
{
	return (a < b) ? b : a;
}

void Camera::setupCameraMatrix()
{
	// Perspective projection matrix
	double projection_matrix[16] = { 0 };
	double aspect_ratio = (double)_camera_film.width / (double)_camera_film.height;
	double fov = pi<double>() / 2.0; // ~ 90 degrees for the field of view
	Transform::perspective(projection_matrix, fov, aspect_ratio, 1.0, 1000.0);

	// Camera view transform
	Transform camera_view;
	camera_view.setPosition(_eye_point.x, _eye_point.y, _eye_point.z);
	camera_view.setRotation(0.0, 90.0, -90.0);

	// Invert the view matrix since it's a camera view
	Transform::invertMatrix(camera_view.matrix, camera_view.matrix);

	// Combine perspective projection and camera view
	double combined_matrix[16] = { 0.0 };
	// Combine the transform matrix for perspective and camera view
	Transform::multiply(projection_matrix, camera_view.matrix, combined_matrix);
	// Invert the transform matrix since we are reversing the process (from normalized device coordinates to world coordinates)
	Transform::invertMatrix(combined_matrix, combined_matrix);

	// Final transform that can be used to convert pixel coordinates to world coordinates
	_transform_matrix = make_mat4(combined_matrix);
}

void Camera::computePixelWidth()
{
	dvec2 p_1 = normalizedPixelCoord(0, 0);
	dvec2 p_2 = normalizedPixelCoord(1, 0);

	_pixel_width = glm::distance(dvec3(_transform_matrix * vec4(p_1.x, p_1.y, 1.0, 1.0)), dvec3(_transform_matrix * vec4(p_2.x, p_2.y, 1.0, 1.0)));
	std::cout << _pixel_width << std::endl;
}

dvec2 Camera::normalizedPixelCoord(const int& x, const int& y)
{
	return dvec2(1.0 - (x / (_camera_film.width  * 0.5)), (y / (_camera_film.height * 0.5)) - 1.0);
}