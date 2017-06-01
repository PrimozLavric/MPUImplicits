#pragma once

#include <glm/glm.hpp>

using namespace glm;

// Structure used to represent a point in the pointcloud
template <typename T>
struct Point {
	tvec3<T> position;
	tvec3<T> normal;

	Point(tvec3<T> position, tvec3<T> normal) {
		this->position = position;
		this->normal = normal;
	}
};