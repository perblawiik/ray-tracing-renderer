#pragma once

#include "glm\glm.hpp"

#include <iostream>
#include <chrono>

#include "Bitmap.h"
#include "TriangleData.h"
#include "Camera.h"
#include "Sphere.h"

using namespace glm;

int main()
{
	const int width = 800;
	const int height = 800;
	const vec3 camera_position(1.0, 0.0, 0.0);
	const vec3 light_position(5.0, 0.0, 4.0);

	TriangleData triangle_objects;
	{
		const double triangleData[] = {
			//Color           Normal		   Vertex 1            Vertex 2           Vertex 3
			//-- FLOOR --//
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   -3.0, 0.0, -2.0,    10.0, 6.0, -2.0,   0.0, 6.0, -2.0,
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   -3.0, 0.0, -2.0,    13.0, 0.0, -2.0,   10.0, 6.0, -2.0,
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   13.0, 0.0, -2.0,    -3.0, 0.0, -2.0,   0.0, -6.0, -2.0,
			1.0, 1.0, 1.0,    0.0, 0.0, 1.0,   13.0, 0.0, -2.0,    0.0, -6.0, -2.0,   10.0, -6.0, -2.0,
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
			//-- FRONT LEFT WALL --//
			0.0, 1.0, 0.5,   -0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   10.0, -6.0, 5.0,   13.0, 0.0, 5.0,
			0.0, 1.0, 0.5,   -0.894427, 0.447214, 0.0,  10.0, -6.0, -5.0,   13.0, 0.0, -5.0,   13.0, 0.0, 5.0,
			//-- FRONT RIGHT WALL --//
			1.0, 1.0, 0.0,   -0.894427, -0.447214, 0.0,  10.0, 6.0, -5.0,   10.0, 6.0, 5.0,    13.0, 0.0, 5.0,
			1.0, 1.0, 0.0,   -0.894427, -0.447214, 0.0,  13.0, 0.0, 5.0,    13.0, 0.0, -5.0,   10.0, 6.0, -5.0,
			//-- BACK LEFT WALL --//
			0.0, 0.0, 1.0,   0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, 5.0,   -3.0, 0.0, 5.0,
			0.0, 0.0, 1.0,   0.894427, 0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, -6.0, -5.0,   0.0, -6.0, 5.0,
			//-- BACK RIGHT WALL --//
			1.0, 0.0, 1.0,   0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,   0.0, 6.0, 5.0,   0.0, 6.0, -5.0,
			1.0, 0.0, 1.0,   0.894427, -0.447214, 0.0,   -3.0, 0.0, -5.0,   -3.0, 0.0, 5.0,   0.0, 6.0, 5.0
		};
		triangle_objects.loadData(triangleData, 20);
	}
	
	Sphere sphere(vec3(1.0, 1.0, 1.0), vec3(5.0, 0.0, -1.0), 1.0);

	Camera camera(width, height, camera_position);
	camera.loadSceneObjects(&triangle_objects, &sphere);
	camera.setLightPosition(light_position);

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
