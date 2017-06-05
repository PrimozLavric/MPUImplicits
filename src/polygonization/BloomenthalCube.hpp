#pragma once
#include <glm/glm.hpp>
#include <vector>
#include "../structures/Point.hpp"

using namespace std;
using namespace glm;

namespace Bloomenthal {

	template <class T, class LambdaF>
	class Cube {
	private:
		// Returns distance from the surface of the polygonized object
		LambdaF distanceFunction;

		// Position of the bottom left corner and edge size of the cube
		tvec3<T> position;
		T size;

		// Cube vertices
		vector<tvec3<T> > vertices;
		vector<T> verticesFval;

	public:
		Cube(tvec3<T> position, T size, vector<T> verticesFval, LambdaF distanceFunction);
		~Cube() {};

		void tetrahedralDecomposition(vector<Point<T> > &points);

	private:
		void tetrahedra(vector<Point<T> > &points, size_t a, size_t b, size_t c, size_t d);

		inline T diff() { return 3 * this->size; }
		tvec3<T> normalid(tvec3<T> position);
		tvec3<T> vertid(size_t a, size_t b);
	};

	template <class T, class LambdaF>
	Cube<T, LambdaF>::Cube(tvec3<T> position, T size, vector<T> verticesFval, LambdaF distanceFunction) : position(position), size(size), verticesFval(verticesFval), distanceFunction(distanceFunction) {
		// Setup cube vertices
		vertices.push_back(tvec3<T>(position.x, position.y, position.z));
		vertices.push_back(tvec3<T>(position.x, position.y, position.z + size));
		vertices.push_back(tvec3<T>(position.x, position.y + size, position.z));
		vertices.push_back(tvec3<T>(position.x, position.y + size, position.z + size));
		vertices.push_back(tvec3<T>(position.x + size, position.y, position.z));
		vertices.push_back(tvec3<T>(position.x + size, position.y, position.z + size));
		vertices.push_back(tvec3<T>(position.x + size, position.y + size, position.z));
		vertices.push_back(tvec3<T>(position.x + size, position.y + size, position.z + size));
	}

	template <class T, class LambdaF>
	void Cube<T, LambdaF>::tetrahedralDecomposition(vector<Point<T> > &points) {
		tetrahedra(points, 0, 2, 4, 1);
		tetrahedra(points, 6, 2, 1, 4);
		tetrahedra(points, 6, 2, 3, 1);
		tetrahedra(points, 6, 4, 1, 5);
		tetrahedra(points, 6, 1, 3, 5);
		tetrahedra(points, 6, 3, 7, 5);
	}

	template <class T, class LambdaF>
	void Cube<T, LambdaF>::tetrahedra(vector<Point<T> > &points, size_t a, size_t b, size_t c, size_t d) {
		int index = 0;
		bool apos = false, bpos = false, cpos = false, dpos = false;

		if (this->verticesFval[a] >= 0) {
			apos = true;
			index += 8;
		}
		if (this->verticesFval[b] >= 0) {
			bpos = true;
			index += 4;
		}
		if (this->verticesFval[c] >= 0) {
			cpos = true;
			index += 2;
		}
		if (this->verticesFval[d] >= 0) {
			dpos = true;
			index += 1;
		}

		tvec3<T> e1, e2, e3, e4, e5, e6;
		tvec3<T> n1, n2, n3, n4, n5, n6;

		if (apos != bpos) {
			e1 = vertid(a, b);
			n1 = normalid(e1);
		}
		if (apos != cpos) {
			e2 = vertid(a, c);
			n2 = normalid(e2);
		}
		if (apos != dpos) {
			e3 = vertid(a, d);
			n3 = normalid(e3);
		}
		if (bpos != cpos) {
			e4 = vertid(b, c);
			n4 = normalid(e4);
		}
		if (bpos != dpos) {
			e5 = vertid(b, d);
			n5 = normalid(e5);
		}
		if (cpos != dpos) {
			e6 = vertid(c, d);
			n6 = normalid(e6);
		}

		switch (index) {
		case 1:
			points.push_back(Point<T>(e5, n5));
			points.push_back(Point<T>(e6, n6));
			points.push_back(Point<T>(e3, n3));
			break;
		case 2:
			points.push_back(Point<T>(e2, n2));
			points.push_back(Point<T>(e6, n6));
			points.push_back(Point<T>(e4, n4));
			break;
		case 3:
			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e5, n5));
			points.push_back(Point<T>(e4, n4));

			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e4, n4));
			points.push_back(Point<T>(e2, n2));
			break;
		case 4:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e4, n4));
			points.push_back(Point<T>(e5, n5));
			break;
		case 5:
			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e4, n4));

			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e4, n4));
			points.push_back(Point<T>(e6, n6));
			break;
		case 6:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e2, n2));
			points.push_back(Point<T>(e6, n6));

			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e6, n6));
			points.push_back(Point<T>(e5, n5));
			break;
		case 7:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e2, n2));
			points.push_back(Point<T>(e3, n3));
			break;
		case 8:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e2, n2));
			break;
		case 9:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e5, n5));
			points.push_back(Point<T>(e6, n6));

			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e6, n6));
			points.push_back(Point<T>(e2, n2));
			break;
		case 10:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e6, n6));

			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e6, n6));
			points.push_back(Point<T>(e4, n4));
			break;
		case 11:
			points.push_back(Point<T>(e1, n1));
			points.push_back(Point<T>(e5, n5));
			points.push_back(Point<T>(e4, n4));
			break;
		case 12:
			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e2, n2));
			points.push_back(Point<T>(e4, n4));

			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e4, n4));
			points.push_back(Point<T>(e5, n5));
			break;
		case 13:
			points.push_back(Point<T>(e6, n6));
			points.push_back(Point<T>(e2, n2));
			points.push_back(Point<T>(e4, n4));
			break;
		case 14:
			points.push_back(Point<T>(e5, n5));
			points.push_back(Point<T>(e3, n3));
			points.push_back(Point<T>(e6, n6));
			break;
		}
	}

	template <class T, class LambdaF>
	tvec3<T> Cube<T, LambdaF>::normalid(tvec3<T> position) {
		T fval = distanceFunction(position);
		tvec3<T> v(0, 0, 0);
		v.x = distanceFunction(tvec3<T>(position.x + diff(), position.y, position.z)) - fval;
		v.y = distanceFunction(tvec3<T>(position.x, position.y + diff(), position.z)) - fval;
		v.z = distanceFunction(tvec3<T>(position.x, position.y, position.z + diff())) - fval;

		return normalize(v);
	}

	template <class T, class LambdaF>
	tvec3<T> Cube<T, LambdaF>::vertid(size_t a, size_t b) {

		tvec3<T> pos, neg;
		if (this->verticesFval[a] < 0) {
			pos = this->vertices[b];
			neg = this->vertices[a];
		}
		else {
			pos = this->vertices[a];
			neg = this->vertices[b];
		}

		tvec3<T> vertex;

		for (int i = 0; i < 5; i++) {
			vertex = (pos + neg) * ((T) 0.5);

			// Check on which side of the surface is the point located
			if (distanceFunction(vertex) > 0.0) {
				pos = vertex;
			}
			else {
				neg = vertex;
			}
		}

		return (pos + neg) * ((T) 0.5);
	}
}