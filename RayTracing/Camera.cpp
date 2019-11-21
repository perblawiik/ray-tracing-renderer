#include "Camera.h"

Camera::Camera(const size_t width, const size_t height, const dvec3& eye)
	: _eye_point(eye), _camera_film(width, height), _scene(nullptr),
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
	const int max_reflections = 7;

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

				// Generate a random offset in the pixel plane (used for ray randomization to reduce aliasing)
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
	IntersectionPoint surface_point;

	// Test the ray against all geometries in the scene (function returns true if valid intersection is found)
	if (geometryIntersectionTest(ray, surface_point) == false) {
		return dvec3(0.0);
	}

	// Light source
	if (surface_point.material->surface_type == Material::SurfaceType::LightSource) {
		return surface_point.material->color;
	}
	// Specular surface
	else if (surface_point.material->surface_type == Material::SurfaceType::Specular) {
		// Compute reflected ray
		dvec3 reflect_direction = reflect(ray.direction, surface_point.normal);
		// Add bias to offset ray origin from surface
		//dvec3 bias = ray.direction * EPSILON;
		Ray reflected_ray(surface_point.position, reflect_direction);

		// Recursive path tracing (perfect reflection)
		return tracePath(reflected_ray, reflection_count - 1);
	}
	// Transparent surface
	else if (surface_point.material->surface_type == Material::SurfaceType::Transparent) {
		// Refractive indices 
		//*** This should probably be apart of the material class ***//
		double n_1 = 1.0; // Air
		double n_2 = 1.5; // Glass

		// Check the current medium the ray is traveling in
		double cos_theta = dot(surface_point.normal, ray.direction);
		dvec3 surface_normal = surface_point.normal;
		
		// If theta < 0.0, the normal is facing the opposite direction from incoming ray (current medium is air)
		if (cos_theta < 0.0) {
			cos_theta = -cos_theta;
		}
		else { // Current medium is glass
			std::swap(n_1, n_2);
			surface_normal = -surface_normal;
		}

		// Compute Fresnel's equation for radiance distribution over reflect and refracted ray
		double reflection_ratio = fresnelsEquation(n_1, n_2, cos_theta);

		// If the current medium is thicker than the outgoing medium and refraction ratio == 1, we have total reflection and no refraction
		dvec3 refracted_light(0.0);
		if (reflection_ratio < 1.0)  {
			// Compute the refracted ray
			double n_1_div_n_2 = n_1 / n_2;
			double A = 1.0 - (n_1_div_n_2 * n_1_div_n_2) * (1.0 - cos_theta * cos_theta);

			dvec3 refracted_direction(
				ray.direction * n_1_div_n_2 +
				surface_normal * (n_1_div_n_2 * cos_theta - sqrt(A))
			);

			// Add bias to offset ray origin from the surface
			dvec3 bias = refracted_direction * EPSILON;
			Ray refracted_ray(surface_point.position + bias, refracted_direction);

			// Recursion (multiplied with the refraction distribution coefficient)
			refracted_light = tracePath(refracted_ray, reflection_count - 1) * (1.0 - reflection_ratio);
		}

		// Compute reflected ray
		dvec3 reflect_direction = reflect(ray.direction, surface_normal);

		// Add bias to offset ray origin from surface
		dvec3 bias = reflect_direction * EPSILON;
		Ray reflected_ray(surface_point.position + bias, reflect_direction);

		// Recursion (multiplied with the reflection distribution coefficient)
		dvec3 reflected_light = tracePath(reflected_ray, reflection_count - 1) * reflection_ratio;

		// Final light
		return refracted_light + reflected_light;
	}
	// Diffuse surface
	else if (surface_point.material->surface_type == Material::SurfaceType::Diffuse) {
		// Let the first random number be equal to cos(theta)
		double cos_theta = _distribution(_generator); 

		// Apply Russian roulette to terminate rays.
		double phi = (cos_theta * 2.0 * PI) / surface_point.material->reflection_coefficient;
		if (phi > 2.0 * PI) {
			return dvec3(0.0);
		}

		// Generate and compute a random direction in the hemisphere (note that cos theta is randomly generated)
		dvec3 sample_direction = hemisphereSampleDirection(cos_theta, surface_point.normal);

		// Create a ray from the sample direction
		//dvec3 bias = sample_direction * EPSILON;
		Ray sample_ray(surface_point.position, sample_direction);

		// Bidirectional Reflectance Distribution Function
		dvec3 brdf = surface_point.material->brdf(surface_point.normal, ray.direction, sample_ray.direction);

		// Recursion for the sample ray
		// First note: the variable phi contains the same parts given by dividing with pdf, multiplying with cos(theta),
		// and also division by the reflection coefficient to compensate for the probability of the russian roulette termination.
		// Second note: the pdf is constant (1.0 / (2.0 * PI)) in this case because all of the random directions share the same probability (equiprobability)
		dvec3 indirect_light = tracePath(sample_ray, reflection_count - 1) * brdf * phi;

		// Compute direct light from light source
		dvec3 direct_light = computeDirectLight(surface_point, brdf, 20);

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
		if (intersection_found && t < light_distance && 
			(tetrahedron->material->surface_type != Material::SurfaceType::Transparent)) {
			return true;
		}
	}

	for (Sphere* sphere : _scene->spheres) {
		double d_near, d_far = -1.0;
		bool intersection_found = sphere->rayIntersection(light_ray, d_near, d_far);

		// Shadow
		if (intersection_found && (d_near < light_distance) && 
			(sphere->material->surface_type != Material::SurfaceType::Transparent)) {
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

			//dvec3 bias = point_to_light_direction * 1e-8;
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

dvec3 Camera::hemisphereSampleDirection(const double &cos_theta, const dvec3& surface_normal)
{
	// Theta is the inclination angle 
	double sin_theta = sqrt(1 - cos_theta * cos_theta);
	// Phi is the azimuth angle
	double phi = 2 * PI * _distribution(_generator); // Second random number

	double x = sin_theta * cos(phi);
	double z = sin_theta * sin(phi);

	// We assume that the first random value is cos(theta), which is equal to the y-coordiante
	dvec3 sample_direction_local(x, cos_theta, z);

	// Create a local coordinate system for the hemisphere of the intersection point with the surface normal as the y-axis
	dvec3 local_x_axis(0.0);
	dvec3 local_z_axis(0.0);
	createLocalCoordinateSystem(surface_normal, local_x_axis, local_z_axis);

	// Transform direction from local to world by multiplying with the local coordinate system matrix
	// Note that we only transform the direction so the translation part is not needed (in that case we would use a 4x4 matrix)
	return dvec3(
		sample_direction_local.x * local_z_axis.x + sample_direction_local.y * surface_normal.x + sample_direction_local.z * local_x_axis.x,
		sample_direction_local.x * local_z_axis.y + sample_direction_local.y * surface_normal.y + sample_direction_local.z * local_x_axis.y,
		sample_direction_local.x * local_z_axis.z + sample_direction_local.y * surface_normal.z + sample_direction_local.z * local_x_axis.z
	);
}

bool Camera::geometryIntersectionTest(const Ray& ray, IntersectionPoint& surface_point)
{
	// Test the ray against all triangle objects
	triangleIntersectionTests(ray, surface_point);

	// Test the ray against the sphere
	sphereIntersectionTests(ray, surface_point);
	
	// If distance is negative, no intersection was found in the ray direction
	return !(surface_point.distance < 0.0);
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

void Camera::sphereIntersectionTests(const Ray& ray, IntersectionPoint& closest_point)
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

void Camera::setupCameraMatrix()
{
	// Perspective projection matrix
	double projection_matrix[16] = { 0 };
	double aspect_ratio = (double)_camera_film.width / (double)_camera_film.height;
	double fov = 75.0 * pi<double>() / 180.0; // 75 degrees for the field of view
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
}

inline dvec3 Camera::barycentricToWorldCoordinates(const Triangle& triangle, const double& u, const double& v)
{
	dvec3 u_vec(u);
	dvec3 v_vec(v);
	dvec3 one_vec(1.0);
	return ((one_vec - u_vec - v_vec) * triangle.vertices[0]) + (u_vec * triangle.vertices[1]) + (v_vec * triangle.vertices[2]);
}

inline dvec2 Camera::normalizedPixelCoord(const int& x, const int& y)
{
	return dvec2(1.0 - (x / (_camera_film.width  * 0.5)), (y / (_camera_film.height * 0.5)) - 1.0);
}

inline double Camera::schlicksEquation(const double& n_1, const double& n_2, const double& cos_theta)
{
	double R_0 = ((n_1 - n_2) / (n_1 + n_2)) * ((n_1 - n_2) / (n_1 + n_2));
	return R_0 + (1.0 - R_0) * pow(1.0 - abs(cos_theta), 5);
}

inline double Camera::fresnelsEquation(const double& n_1, const double& n_2, const double& cos_theta1)
{
	// Compute sin_theta2 using Snell's law
	double sin_theta2 = n_1 / n_2 * sqrt(max(0.0, 1.0 - cos_theta1 * cos_theta1));

	// If larger or equal to one, we have total reflection (internal reflection inside a thicker medium compared to outside)
	if (sin_theta2 >= 1.0) {
		return 1.0;
	}

	// Compute cos_theta2 (derived from cos_theta^2 + sin_theta^2 = 1)
	double cos_theta2 = sqrt(max(0.0, 1.0 - sin_theta2 * sin_theta2));

	// Calculate the square roots of the Fresnel equations
	double Rs_sqroot = ((n_1 * cos_theta1) - (n_2 * cos_theta2)) / ((n_1 * cos_theta1) + (n_2 * cos_theta2));
	double Rp_sqroot = ((n_2 * cos_theta1) - (n_1 * cos_theta2)) / ((n_2 * cos_theta1) + (n_1 * cos_theta2));

	// The reflection ratio is given by the average of the two equations
	return (Rs_sqroot * Rs_sqroot + Rp_sqroot * Rp_sqroot) * 0.5;
}