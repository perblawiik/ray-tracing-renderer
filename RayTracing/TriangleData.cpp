#include "TriangleData.h"

TriangleData::TriangleData() {}

TriangleData::TriangleData(const double data[], const size_t num_triangles)
{
	loadData(data, num_triangles);
}

void TriangleData::loadData(const double data[], const size_t num_triangles)
{
	// Clear previous data
	triangles.clear();
	// Reserve capacity for the triangles
	triangles.reserve(num_triangles);
	// Triangle data per row
	size_t stride = 15;
	for (size_t i = 0; i < num_triangles; ++i) {
		triangles.emplace_back(
			Triangle(
				glm::vec3(data[i * stride], data[i * stride + 1], data[i * stride + 2]),  // Color
				glm::vec3(data[i * stride + 3], data[i * stride + 4], data[i * stride + 5]),  // Normal
				glm::vec3(data[i * stride + 6], data[i * stride + 7], data[i * stride + 8]),  // Vertex 1
				glm::vec3(data[i * stride + 9], data[i * stride + 10], data[i * stride + 11]), // Vertex 2
				glm::vec3(data[i * stride + 12], data[i * stride + 13], data[i * stride + 14])  // Vertex 3
			)
		);
	}
}

