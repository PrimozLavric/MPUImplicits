#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "../structures/Point.hpp"

using namespace std;
using namespace glm;

namespace MPUIUtility {
	template <typename T>
	inline T weight(tvec3<T> position, tvec3<T> center, T radius) {
		position -= center;
		// Distance between the point and the center
		T dist = sqrt(dot(position, position));

		if (dist > radius) {
			return 0;
		}

		//B-spline (degree = 2)
		dist = 1.5 * (dist / radius);
		if (dist < 0.5) {
			return (-dist * dist + 0.75);
		}
		else {
			return (0.5 * (1.5 - dist) * (1.5 - dist));
		}
	}

	template <typename T>
	inline vector<tvec3<T> > generateAuxilaryPoints(tvec3<T> center, T cellSize) {
		vector<tvec3<T> > auxP;

		T hSize = cellSize / 2;
		auxP.push_back(tvec3<T>(center.x - hSize, center.y - hSize, center.z - hSize));
		auxP.push_back(tvec3<T>(center.x + hSize, center.y - hSize, center.z - hSize));
		auxP.push_back(tvec3<T>(center.x - hSize, center.y + hSize, center.z - hSize));
		auxP.push_back(tvec3<T>(center.x + hSize, center.y + hSize, center.z - hSize));
		auxP.push_back(tvec3<T>(center.x - hSize, center.y - hSize, center.z + hSize));
		auxP.push_back(tvec3<T>(center.x + hSize, center.y - hSize, center.z + hSize));
		auxP.push_back(tvec3<T>(center.x - hSize, center.y + hSize, center.z + hSize));
		auxP.push_back(tvec3<T>(center.x + hSize, center.y + hSize, center.z + hSize));
		auxP.push_back(tvec3<T>(center.x, center.y, center.z));

		return auxP;
	}


	template <typename T>
	tvec3<T> positionAverage(vector<Point<T> > points) {
		tvec3<T> avgPoint(0, 0, 0);

		for (auto it = points.begin(); it != points.end(); it++) {
			avgPoint += it->position;
		}

		return avgPoint / ((T) points.size());
	}

	template <typename T>
	tvec3<T> normalAverage(vector<Point<T> > points) {
		tvec3<T> avgNormal(0, 0, 0);

		for (auto it = points.begin(); it != points.end(); it++) {
			avgNormal += it->normal;
		}

		return avgNormal / ((T) points.size());
	}
}
