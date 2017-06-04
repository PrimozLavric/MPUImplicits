#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "../structures/Point.hpp"

using namespace std;
using namespace glm;

namespace MPUIUtility {
	template <typename T>
	inline T weight(tvec3<T> position, tvec3<T> center, T radius) {
		tvec3<T> pc = position - center;
		// Distance between the point and the center
		T dist = sqrt(dot(pc, pc));

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
	pair<tvec3<T>, tvec3<T> > findAveragedTangentPlane(vector<Point<T> > points, tvec3<T> cellCenter, T radius) {
		// First - position, second - normal
		pair<tvec3<T>, tvec3<T> > tangentP = make_pair(tvec3<T>(0, 0, 0), tvec3<T>(0, 0, 0));

		T weightSum = 0;
		for (auto itP = points.begin(); itP != points.end(); itP++) {
			T w = weight(itP->position, cellCenter, radius);
			weightSum += w;

			tangentP.first += itP->position * w;
			tangentP.second += itP->normal * w;
		}

		tangentP.first /= weightSum;
		
		// Normalize the normal
		if (length(tangentP.second) == 0) {
			tangentP.second = tvec3<T>(1, 1, 1);
		}
		
		tangentP.second = normalize(tangentP.second);

		return tangentP;
	}
}
