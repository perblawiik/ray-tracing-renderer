#pragma once

#include "Material.h"

class Object
{
public:
	Object(const Material& material_);
	bool intersectionTest()
	{
	};

	Material material;
};