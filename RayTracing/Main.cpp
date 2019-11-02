#pragma once

#include "glm\glm.hpp"

#include <iostream>
#include <chrono>
#include <string>

#include "Bitmap.h"
#include "TriangleObject.h"
#include "Camera.h"
#include "Sphere.h"
#include "AreaLightSource.h"
#include "Material.h"
#include "Scene.h"

using namespace glm;

int main()
{
	const int width = 800;
	const int height = 800;
	const dvec3 camera_position(12.0, 0.0, 0.0);

	Scene scene;

	AreaLightSource area_light_source(
		dvec3(0.0, 1.0, 4.999995), // V0
		dvec3(2.0, 0.0, 4.999995), // V1
		dvec3(0.0, -1.0, 4.999995), // V2
		dvec3(1.0), // Color
		100.0 // Watts
	);
	scene.addLightSource(&area_light_source);

	TriangleObject floor(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0), 0.75));
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

	TriangleObject ceiling(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0), 0.75));
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

	TriangleObject left_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0, 0.25, 0.25), 0.75));
	{
		const double triangle_data[] = {
			// Normal		  Vertex 1           Vertex 2           Vertex 3
			0.0, 1.0, 0.0,    0.0, -6.0, -5.0,   10.0, -6.0, 5.0,    0.0, -6.0, 5.0,
			0.0, 1.0, 0.0,    0.0, -6.0, -5.0,   10.0, -6.0, -5.0,  10.0, -6.0, 5.0
		};
		left_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&left_wall);

	TriangleObject right_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(0.0, 1.0, 1.0), 0.75));
	{
		const double triangle_data[] = {
			// Normal		  Vertex 1           Vertex 2           Vertex 3
			0.0, -1.0, 0.0,   0.0, 6.0, -5.0,    0.0, 6.0, 5.0,     10.0, 6.0, 5.0,
			0.0, -1.0, 0.0,   0.0, 6.0, -5.0,    10.0, 6.0, 5.0,    10.0, 6.0, -5.0
		};
		right_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&right_wall);

	TriangleObject back_right_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(0.0, 1.0, 0.5), 0.75));
	{
		const double triangle_data[] = {
			// Normal		           Vertex 1            Vertex 2           Vertex 3
			-0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   10.0, -6.0, 5.0,   13.0, 0.0, 5.0,
			-0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   13.0, 0.0, -5.0,   13.0, 0.0, 5.0
		};
		back_right_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&back_right_wall);

	TriangleObject back_left_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0, 0.0, 1.0), 0.75));
	{
		const double triangle_data[] = {
			// Normal		            Vertex 1           Vertex 2           Vertex 3
			-0.894427, -0.447214, 0.0,  10.0, 6.0, -5.0,   10.0, 6.0, 5.0,    13.0, 0.0, 5.0,
			-0.894427, -0.447214, 0.0,  13.0, 0.0, 5.0,    13.0, 0.0, -5.0,   10.0, 6.0, -5.0
		};
		back_left_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&back_left_wall);

	TriangleObject front_left_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(0.0, 0.0, 1.0), 0.75));
	{
		const double triangle_data[] = {
			// Normal		           Vertex 1           Vertex 2           Vertex 3
			0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, 5.0,   -3.0, 0.0, 5.0,
			0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, -5.0,   0.0, -6.0, 5.0
		};
		front_left_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&front_left_wall);

	TriangleObject front_right_wall(new Lambertian(Material::SurfaceType::Diffuse, dvec3(1.0, 1.0, 0.0), 0.75));
	{
		const double triangle_data[] = {
			// Normal		            Vertex 1           Vertex 2         Vertex 3
			0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, 6.0, 5.0,   0.0, 6.0, -5.0,
			0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,  -3.0, 0.0, 5.0,   0.0, 6.0, 5.0
		};
		front_right_wall.loadData(triangle_data, 2);
	}
	scene.addTriangleObject(&front_right_wall);

	TriangleObject specular_tetrahedron(new Lambertian(Material::SurfaceType::Specular, dvec3(1.0, 1.0, 1.0), 1.0));
	{
		dvec3 position(0.0, -4.0, -5.0);
		double scale = 4.0;
		specular_tetrahedron.createTetrahedron(position, scale);
	}
	scene.addTriangleObject(&specular_tetrahedron);
	
	Lambertian* specular_sphere_material = new Lambertian(Material::SurfaceType::Specular, dvec3(1.0, 1.0, 1.0), 1.0);
	Sphere specular_sphere(specular_sphere_material, dvec3(3.0, 1.0, -2.5), 2.0);
	scene.addSphere(&specular_sphere);

	Camera camera(width, height, camera_position);
	camera.loadScene(&scene);

	auto time_start = std::chrono::high_resolution_clock::now();

	std::cout << "Rendering..." << std::endl;

	// Render camera view
	const int num_samples = 32;
	camera.render(num_samples);

	auto time_end = std::chrono::high_resolution_clock::now();
	auto run_time = std::chrono::duration<double, std::milli>(time_end - time_start).count();

	std::string file_name = "render_N" + std::to_string(num_samples) + "_" + std::to_string((int)(run_time / 1000.0)) + "s.bmp";

	// Create and save image file
	camera.createImage(file_name.c_str());

	//std::cout << "Render complete: " << run_time / 1000.0 << " seconds" << std::endl;

	//system("pause");

	return 0;
}
