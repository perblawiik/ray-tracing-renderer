#include "TriangleObject.h"

TriangleObject::TriangleObject(Material* material) 
	: material(material) 
{}

void TriangleObject::loadData(const double data[], const size_t num_triangles)
{
	// Clear previous data
	triangles.clear();
	// Reserve capacity for the triangles
	triangles.reserve(num_triangles);
	// Triangle data per row
	size_t stride = 12;
	for (size_t i = 0; i < num_triangles; ++i) {
		triangles.emplace_back(
			Triangle(
				glm::dvec3(data[i * stride], data[i * stride + 1], data[i * stride + 2]),  // Normal
				glm::dvec3(data[i * stride + 3], data[i * stride + 4], data[i * stride + 5]),  // Vertex 1
				glm::dvec3(data[i * stride + 6], data[i * stride + 7], data[i * stride + 8]),  // Vertex 2
				glm::dvec3(data[i * stride + 9], data[i * stride + 10], data[i * stride + 11]) // Vertex 3
			)
		);
	}
}

void TriangleObject::createTetrahedron(const glm::dvec3& pos, const double& scale = 1.0)
{
	const double triangle_data[] = {
		// Normal					Vertex 1                      Vertex 2					    Vertex 3
		// Triangle 1
		0.0, 0.0, -1.0,             pos.x, pos.y, pos.z,		  pos.x, pos.y + scale, pos.z,  pos.x + scale, pos.y, pos.z,
		// Triangle 2
		0.0, -1.0, 0.0,             pos.x, pos.y, pos.z,		  pos.x + scale, pos.y, pos.z,  pos.x, pos.y, pos.z + scale,
		// Triangle 3
		-1.0, 0.0, 1.0,             pos.x, pos.y, pos.z,          pos.x, pos.y, pos.z + scale,  pos.x, pos.y + scale, pos.z,
		// Triangle 4
		0.57735, 0.57735, 0.57735,  pos.x + scale, pos.y, pos.z,  pos.x, pos.y + scale, pos.z,  pos.x, pos.y, pos.z + scale
	};
	loadData(triangle_data, 4);
}

