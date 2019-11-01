#pragma once

#include "glm\glm.hpp"

#include <iostream>
#include <chrono>

#include "Bitmap.h"
#include "TriangleObject.h"
#include "Camera.h"
#include "Sphere.h"
#include "AreaLightSource.h"

using namespace glm;

int main()
{
	const int width = 800;
	const int height = 800;
	const vec3 camera_position(12.0, 0.0, 0.0);
	const vec3 light_position(0.0, 0.0, 4.0);

	AreaLightSource area_light_source(
		vec3(0.0, 1.0, 4.999995), // V0
		vec3(2.0, 0.0, 4.999995), // V1
		vec3(0.0, -1.0, 4.999995), // V2
		dvec3(1.0), // Color
		100.0 // Intensity
	);

	TriangleObject walls;
	{
		const double triangle_data[] = {
			//Color           Normal		   Vertex 1            Vertex 2           Vertex 3
			//-- FLOOR --//
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   -3.0, 0.0, -5.0,    13.0, 0.0, -5.0,   10.0, 6.0, -5.0,
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   -3.0, 0.0, -5.0,    10.0, 6.0, -5.0,   0.0, 6.0, -5.0,
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   13.0, 0.0, -5.0,    -3.0, 0.0, -5.0,   0.0, -6.0, -5.0,
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   13.0, 0.0, -5.0,    0.0, -6.0, -5.0,   10.0, -6.0, -5.0,
			//-- CEILING --//
			1.0, 1.0, 1.0,    0.0, 0.0, -1.0,   -3.0, 0.0, 5.0,    0.0, 6.0, 5.0,     10.0, 6.0, 5.0,
			1.0, 1.0, 1.0,    0.0, 0.0, -1.0,   -3.0, 0.0, 5.0,    10.0, 6.0, 5.0,    13.0, 0.0, 5.0,
			1.0, 1.0, 1.0,    0.0, 0.0, -1.0,   13.0, 0.0, 5.0,    0.0, -6.0, 5.0,    -3.0, 0.0, 5.0,
			1.0, 1.0, 1.0,    0.0, 0.0, -1.0,   13.0, 0.0, 5.0,    10.0, -6.0, 5.0,    0.0, -6.0, 5.0,
			//-- LEFT WALL --//
			1.0, 0.25, 0.25,  0.0, 1.0, 0.0,    0.0, -6.0, -5.0,   10.0, -6.0, 5.0,    0.0, -6.0, 5.0,
			1.0, 0.25, 0.25,  0.0, 1.0, 0.0,    0.0, -6.0, -5.0,   10.0, -6.0, -5.0,  10.0, -6.0, 5.0,
			//-- RIGHT WALL --//
			0.0, 1.0, 1.0,    0.0, -1.0, 0.0,   0.0, 6.0, -5.0,     0.0, 6.0, 5.0,     10.0, 6.0, 5.0,
			0.0, 1.0, 1.0,    0.0, -1.0, 0.0,   0.0, 6.0, -5.0,     10.0, 6.0, 5.0,    10.0, 6.0, -5.0,
			//-- BACK RIGHT WALL --//
			0.0, 1.0, 0.5,   -0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   10.0, -6.0, 5.0,   13.0, 0.0, 5.0,
			0.0, 1.0, 0.5,   -0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   13.0, 0.0, -5.0,   13.0, 0.0, 5.0,
			//-- BACK LEFT WALL --//
			1.0, 0.0, 1.0,   -0.894427, -0.447214, 0.0,  10.0, 6.0, -5.0,   10.0, 6.0, 5.0,    13.0, 0.0, 5.0,
			1.0, 0.0, 1.0,   -0.894427, -0.447214, 0.0,  13.0, 0.0, 5.0,    13.0, 0.0, -5.0,   10.0, 6.0, -5.0,
			//-- FRONT LEFT WALL --//
			0.0, 0.0, 1.0,   0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, 5.0,   -3.0, 0.0, 5.0,
			0.0, 0.0, 1.0,   0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, -5.0,   0.0, -6.0, 5.0,
			//-- FRONT RIGHT WALL --//
			1.0, 1.0, 0.0,   0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, 6.0, 5.0,   0.0, 6.0, -5.0,
			1.0, 1.0, 0.0,   0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,   -3.0, 0.0, 5.0,   0.0, 6.0, 5.0
		};
		walls.loadData(triangle_data, 20);
	}

	TriangleObject tetrahedron;
	{
		dvec3 color(1.0, 1.0, 1.0);
		vec3 position(0.0, -4.0, -5.0);
		double scale = 4.0;
		tetrahedron.createTetrahedron(color, position, scale);
	}
	
	Sphere sphere(vec3(1.0, 1.0, 1.0), vec3(3.0, 1.0, -2.5), 2.0);

	Camera camera(width, height, camera_position);
	camera.loadSceneObjects(&walls, &sphere, &tetrahedron);
	camera.setLightPosition(light_position);
	camera.addLightSource(&area_light_source);

	auto time_start = std::chrono::high_resolution_clock::now();

	std::cout << "Rendering..." << std::endl;

	// Render from camera to bitmap image
	camera.render();

	// Create and save image file
	camera.createImage("render.bmp");

	auto time_end = std::chrono::high_resolution_clock::now();
	auto run_time = std::chrono::duration<double, std::milli>(time_end - time_start).count();

	std::cout << "Render complete: " << run_time / 1000.0 << " seconds" << std::endl;

	system("pause");

	return 0;
}
