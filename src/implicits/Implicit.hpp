#pragma once
#include <glm/glm.hpp>

using namespace glm;

template <class T>
class Implicit {
public:
	virtual T funValue(tvec3<T> position) = 0;

	virtual tvec3<T> gradient(tvec3<T> position) = 0;

	virtual T getAproximationError() = 0;
};