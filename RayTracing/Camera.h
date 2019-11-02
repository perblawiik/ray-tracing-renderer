#pragma once

#include "glm/vec3.hpp"
#include "glm/gtc/constants.hpp"

#include <omp.h>
#include <vector>
#include <random>

#include "TriangleObject.h"
#include "Bitmap.h"
#include "Sphere.h"
#include "AreaLightSource.h"

using namespace glm;

#define PI 3.14159265358979323846

struct IntersectionPoint
{
	enum SurfaceType {
		Diffuse,
		Specular,
		Transparent,
		LightSource
	};

	double distance;
	dvec3 color;
	vec3 position;
	vec3 normal;
	SurfaceType type;

	IntersectionPoint() : distance(-1.0), color(dvec3(0.0, 0.0, 0.0)), position(vec3(0.0, 0.0, 0.0)), normal(vec3(0.0, 0.0, 0.0)), type(SurfaceType::Diffuse){}
};

struct Film
{
	size_t width;
	size_t height;
	vec3 position;
	std::vector<dvec3> pixel_data;

	Film(const size_t& width, const size_t& height, const vec3& position)
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
	Camera(const size_t width, const size_t height, const vec3& eye);

	void setLightPosition(const vec3& position);
	void addLightSource(AreaLightSource* light_source);

	// Saves pointer to all scene objects used when rendering
	void loadSceneObjects(TriangleObject* walls, Sphere* sphere, TriangleObject* tetrahedron);
	
	// Render scene with ray tracing
	void render();

	void createImage(const char* file_name);

private:
	vec3 _eye_point;
	Film _camera_film;

	TriangleObject* _diffuse_walls;
	TriangleObject* _specular_tetrahedron;
	Sphere* _specular_sphere;
	vec3 _light_position;
	AreaLightSource* _light_source;

	std::default_random_engine _generator;
	std::uniform_real_distribution<double> _distribution;

	dvec3 tracePath(const Ray& ray, const int reflection_count);

	bool shadowRay(const vec3& surface_point, const vec3& point_to_light, const double& light_distance);

	dvec3 computeDirectLight(const IntersectionPoint& surface_point, const dvec3& brdf, const size_t sample_ray_count);

	vec3 hemisphereSampleDirection(const double &random_1, const double &random_2, const vec3& surface_normal);

	void createLocalCoordinateSystem(const vec3& N, vec3& Nt, vec3& Nb);

	void triangleIntersectionTests(const Ray& ray, IntersectionPoint& closest_point);

	void sphereIntersectionTest(const Ray& ray, IntersectionPoint& closest_point);

	inline vec3 barycentricToWorldCoordinates(const Triangle& triangle, const double& u, const double& v);

	inline double max(const double& a, const double& b);
};