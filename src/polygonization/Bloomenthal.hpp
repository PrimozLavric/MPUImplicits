#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "../structures/Point.hpp"
#include "BloomenthalCube.hpp"

using namespace std;
using namespace glm;

namespace Bloomenthal {
	
	template <typename T, typename LambdaF>
	vector<Point<T> > polygonize(T size, LambdaF distanceFunction) {

		tvec3<T> cubePos(-0.5, -0.5, -0.5);
		T step = 1.0 / size;

		vector<vector<vector<double>>> verticesFval(size + 1, vector<vector<T>>(size + 1, vector<T>(size + 1)));
		vector<vector<vector<double>>> active(size + 1, vector<vector<T>>(size + 1, vector<T>(size + 1, false)));

		// Calculate distance values for vertices of each cell
		/*for (int i = 0; i <= size; i++) {
			verticesFval.push_back(vector<vector<T> >());
			mask.push_back(vector<vector<bool> >());

			for (int j = 0; j <= size; j++) {
				verticesFval[i].push_back(vector<T>());

				for (int k = 0; k <= size; k++) {
					verticesFval[i][j].push_back(distanceFunction(tvec3<T>(cubePos.x + i * step, cubePos.y + j * step, cubePos.z + k * step)));
				}
			}
		}*/

		// Holds the generated points
		vector<Point<T> > points;

		vector<T> cubeVtxFval(8);
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				for (int k = 0; k < size; k++) {
					// Fetch vertex values for this cube
					cubeVtxFval[0] = verticesFval[i][j][k];
					cubeVtxFval[1] = verticesFval[i][j][k + 1];
					cubeVtxFval[2] = verticesFval[i][j + 1][k];
					cubeVtxFval[3] = verticesFval[i][j + 1][k + 1];
					cubeVtxFval[4] = verticesFval[i + 1][j][k];
					cubeVtxFval[5] = verticesFval[i + 1][j][k + 1];
					cubeVtxFval[6] = verticesFval[i + 1][j + 1][k];
					cubeVtxFval[7] = verticesFval[i + 1][j + 1][k + 1];

					// Generate cube and decompose
					Cube<T, LambdaF> c(tvec3<T>(cubePos.x + i * step, cubePos.y + j * step, cubePos.z + k * step), step, cubeVtxFval, distanceFunction);
					c.tetrahedralDecomposition(points);
				}
			}
		}

		return points;
	}

}