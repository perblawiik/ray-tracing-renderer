#pragma once

#include "glm/vec3.hpp"
#include "glm/mat4x4.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/gtc/constants.hpp"

#include <omp.h>
#include <vector>
#include <random>

#include "Material.h"
#include "TriangleObject.h"
#include "Bitmap.h"
#include "Sphere.h"
#include "AreaLightSource.h"
#include "Scene.h"
#include "Transform.h"

using namespace glm;

#define PI 3.14159265358979323846

struct IntersectionPoint
{
	double distance;
	dvec3 position;
	dvec3 normal;
	Material* material;

	IntersectionPoint() : distance(-1.0), position(dvec3(0.0, 0.0, 0.0)), normal(dvec3(0.0, 0.0, 0.0)), material(nullptr) { }
};

struct Film
{
	size_t width;
	size_t height;
	dvec3 position;
	std::vector<dvec3> pixel_data;

	Film(const size_t& width, const size_t& height, const dvec3& position)
		: width(width), height(height), position(position)
	{
		pixel_data.resize(width * height);
	}

	void addPixelData(const size_t& x, const size_t& y, const dvec3& color) 
	{
		size_t index = y * width + x;
		pixel_data[index] += color;
	}
};

class Camera
{
public:
	// Constructor
	Camera(const size_t width, const size_t height, const dvec3& eye);

	void loadScene(Scene* scene);
	
	// Render scene with ray tracing
	void render(const int& num_samples);

	void createImage(const char* file_name);

private:
	dvec3 _eye_point;
	Film _camera_film;

	Scene* _scene;

	dmat4 _transform_matrix;

	double _pixel_width;

	std::default_random_engine _generator;
	std::uniform_real_distribution<double> _distribution;

	dvec3 tracePath(const Ray& ray, const int reflection_count);

	bool shadowRay(const dvec3& surface_point, const dvec3& point_to_light, const double& light_distance);

	dvec3 computeDirectLight(const IntersectionPoint& surface_point, const dvec3& brdf, const size_t sample_ray_count);

	dvec3 hemisphereSampleDirection(const double &random_1, const double &random_2, const dvec3& surface_normal);

	void createLocalCoordinateSystem(const dvec3& N, dvec3& Nt, dvec3& Nb);

	void triangleIntersectionTests(const Ray& ray, IntersectionPoint& closest_point);

	void sphereIntersectionTest(const Ray& ray, IntersectionPoint& closest_point);

	inline dvec3 barycentricToWorldCoordinates(const Triangle& triangle, const double& u, const double& v);

	inline double max(const double& a, const double& b);

	void setupCameraMatrix();

	void computePixelWidth();

	inline dvec2 normalizedPixelCoord(const int& x, const int& y);

	inline double schlicksEquation(const double& n_1, const double& n_2, const double& cos_theta);
};