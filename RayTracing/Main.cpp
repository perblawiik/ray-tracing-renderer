#pragma once

#include "glm\vec3.hpp"

#include <iostream>
#include <chrono>
#include <string>

#include "Bitmap.h"
#include "TriangleObject.h"
#include "Camera.h"
#include "Sphere.h"
#include "AreaLightSource.h"
#include "Lambertian.h"
#include "Scene.h"

using namespace glm;

const char* createFileName(const std::string& scene_name, const int& width, const int& height, const int& num_samples, const double& run_time);

int main()
{
	//******** 1. Set up the camera ********//
	//--------------------------------------//

	// Pixel dimensions for the image to render
	const int width = 1024;
	const int height = 768;

	const dvec3 camera_position(12.0, 0.0, 0.0);
	Camera camera(width, height, camera_position);

	//******** 2. Set up the scene ********//
	//-------------------------------------//

	// Initiate the scene component that will store all the scene objects
	Scene scene;

	// Create two triangular area light sources to form one rectangular light
	AreaLightSource area_light_source_1(
		dvec3(3.0, 1.0, 4.999995), // V0
		dvec3(5.0, 1.0, 4.999995), // V1
		dvec3(3.0, -1.0, 4.999995), // V2
		new Lambertian(Material::SurfaceType::LightSource, dvec3(1.0), 1.0), // Material
		50 // Watts
	);
	scene.addLightSource(&area_light_source_1);

	AreaLightSource area_light_source_2(
		dvec3(5.0, 1.0, 4.999995), // V0
		dvec3(5.0, -1.0, 4.999995), // V1
		dvec3(3.0, -1.0, 4.999995), // V2
		new Lambertian(Material::SurfaceType::LightSource, dvec3(1.0), 1.0), // Material
		50.0 // Watts
	);
	scene.addLightSource(&area_light_source_2);
	
	// Define the reflection coefficient used on all diffuse surfaces (walls, ceiling, floor) 
	double diffuse_reflection_coefficient = 0.75;

	// Define the room (walls, ceiling, floor)
	TriangleObject floor(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		 Vertex 1            Vertex 2           Vertex 3
			0.0, 0.0, 1.0,   -3.0, 0.0, -5.0,    13.0, 0.0, -5.0,   10.0, 6.0, -5.0,
			0.0, 0.0, 1.0,   -3.0, 0.0, -5.0,    10.0, 6.0, -5.0,   0.0, 6.0, -5.0,
			0.0, 0.0, 1.0,   13.0, 0.0, -5.0,    -3.0, 0.0, -5.0,   0.0, -6.0, -5.0,
			0.0, 0.0, 1.0,   13.0, 0.0, -5.0,    0.0, -6.0, -5.0,   10.0, -6.0, -5.0
		};
		floor.loadData(triangle_data, 4);
	}
	scene.addTriangleObject(&floor);

	TriangleObject ceiling(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		  Vertex 1           Vertex 2           Vertex 3
			0.0, 0.0, -1.0,   -3.0, 0.0, 5.0,    0.0, 6.0, 5.0,     10.0, 6.0, 5.0,
			0.0, 0.0, -1.0,   -3.0, 0.0, 5.0,    10.0, 6.0, 5.0,    13.0, 0.0, 5.0,
			0.0, 0.0, -1.0,   13.0, 0.0, 5.0,    0.0, -6.0, 5.0,    -3.0, 0.0, 5.0,
			0.0, 0.0, -1.0,   13.0, 0.0, 5.0,    10.0, -6.0, 5.0,    0.0, -6.0, 5.0
		};
		ceiling.loadData(triangle_data, 4);
	}
	scene.addTriangleObject(&ceiling);

	TriangleObject left_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0, 0.25, 0.25), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		  Vertex 1           Vertex 2           Vertex 3
			0.0, 1.0, 0.0,    0.0, -6.0, -5.0,   10.0, -6.0, 5.0,    0.0, -6.0, 5.0,
			0.0, 1.0, 0.0,    0.0, -6.0, -5.0,   10.0, -6.0, -5.0,  10.0, -6.0, 5.0
		};
		left_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&left_wall);

	TriangleObject right_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(0.0, 1.0, 1.0), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		  Vertex 1           Vertex 2           Vertex 3
			0.0, -1.0, 0.0,   0.0, 6.0, -5.0,    0.0, 6.0, 5.0,     10.0, 6.0, 5.0,
			0.0, -1.0, 0.0,   0.0, 6.0, -5.0,    10.0, 6.0, 5.0,    10.0, 6.0, -5.0
		};
		right_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&right_wall);

	TriangleObject back_right_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(0.0, 1.0, 0.5), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		           Vertex 1            Vertex 2           Vertex 3
			-0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   10.0, -6.0, 5.0,   13.0, 0.0, 5.0,
			-0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   13.0, 0.0, -5.0,   13.0, 0.0, 5.0
		};
		back_right_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&back_right_wall);

	TriangleObject back_left_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0, 0.0, 1.0), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		            Vertex 1           Vertex 2           Vertex 3
			-0.894427, -0.447214, 0.0,  10.0, 6.0, -5.0,   10.0, 6.0, 5.0,    13.0, 0.0, 5.0,
			-0.894427, -0.447214, 0.0,  13.0, 0.0, 5.0,    13.0, 0.0, -5.0,   10.0, 6.0, -5.0
		};
		back_left_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&back_left_wall);

	TriangleObject front_left_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(0.157, 0.208, 0.576), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		           Vertex 1           Vertex 2           Vertex 3
			0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, 5.0,   -3.0, 0.0, 5.0,
			0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, -5.0,   0.0, -6.0, 5.0
		};
		front_left_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&front_left_wall);

	TriangleObject front_right_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0, 1.0, 0.0), diffuse_reflection_coefficient));
	{
		const double triangle_data[] = {
			// Normal		            Vertex 1           Vertex 2         Vertex 3
			0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, 6.0, 5.0,   0.0, 6.0, -5.0,
			0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,  -3.0, 0.0, 5.0,   0.0, 6.0, 5.0
		};
		front_right_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&front_right_wall);

	// Define a specular tetrahedron
	TriangleObject specular_tetrahedron(new Lambertian(Material::SurfaceType::Specular, dvec3(1.0, 1.0, 1.0), 1.0));
	{
		dvec3 position(0.0, -6.0, -5.0);
		double scale = 4.0;
		specular_tetrahedron.createTetrahedron(position, scale);
	}
	scene.addTriangleObject(&specular_tetrahedron);
	
	// Define a specular sphere
	Lambertian* specular_sphere_material = new Lambertian(Material::SurfaceType::Specular, dvec3(1.0, 1.0, 1.0), 1.0);
	Sphere specular_sphere(specular_sphere_material, dvec3(2.5, 4.0, -2.0), 2.0);
	scene.addSphere(&specular_sphere);

	// Define a transparent sphere
	Lambertian* transparent_sphere_material = new Lambertian(Material::SurfaceType::Transparent, dvec3(1.0, 1.0, 1.0), 1.0);
	Sphere transparent_sphere(transparent_sphere_material, dvec3(3.25, -0.5, -3.0), 2.0);
	scene.addSphere(&transparent_sphere);
	
	// Add a pointer to the scene component to the camera
	camera.loadScene(&scene);

	//******** 3. Render the scene ********//
	//-------------------------------------//

	auto time_start = std::chrono::high_resolution_clock::now();

	std::cout << "Rendering..." << std::endl;

	// Render camera view
	const int num_samples = 512;
	camera.render(num_samples);

	auto time_end = std::chrono::high_resolution_clock::now();
	auto run_time = std::chrono::duration<double, std::milli>(time_end - time_start).count();

	const char* file_name = createFileName("scene4", width, height, num_samples, run_time);
	
	// Create and save image file
	camera.createImage(file_name);

	//std::cout << "Render complete: " << run_time / 1000.0 << " seconds" << std::endl;
	//system("pause");

	return 0;
}

const char* createFileName(const std::string& scene_name, const int& width, const int& height, const int& num_samples, const double& run_time)
{
	std::string resolution = std::to_string(width) + "x" + std::to_string(height);
	std::string samples = "N" + std::to_string(num_samples);
	std::string time = std::to_string((int)(run_time / 1000.0)) + "s";
	std::string file_name = "Render/" + scene_name + "_" + resolution + "_" + samples + "_" + time + ".bmp";
	return file_name.c_str();
}
