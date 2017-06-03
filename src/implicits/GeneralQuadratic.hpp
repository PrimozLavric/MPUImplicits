#pragma once
#include "Implicit.hpp"
#include <Eigen/Dense>
#include <glm/glm.hpp>
#include "../structures/Point.hpp"
#include "UtilityFunctions.hpp"
#include <vector>
#include <algorithm>

using namespace std;
using namespace glm;
using namespace Eigen;


template <class T>
class GeneralQuadratic : public Implicit<T> {
private:
	// Holds the coefficients
	T c[10];

	// True if the coefficients are initialized
	bool initialized;

	// Number of nearest neighbours for auxilary points
	const size_t auxN;

	// TODO: Perhaps it would be better if the weight was set based on the number of accepted auxilary points
	const T auxWeight = 1.0 / 9.0;

	// Threshold for when the singular value is considered 0
	const T svThresh;

	T maxError;

public:
	// Fits the quadric implicit based on the support points, auxilary points and cell center
	GeneralQuadratic(size_t auxN = 6, T svThresh = 1E-7) : initialized(false), auxN(auxN), svThresh(svThresh) {};
	~GeneralQuadratic() {};

	bool fitOnPoints(vector<Point<T>> &points, tvec3<T> cellCenter, T cellSize, T radius);

	inline tvec3<T> gradient(tvec3<T> position);

	inline T funValue(tvec3<T> position);

	inline T getAproximationError() { return this->maxError; }

private:
	void calculateError(vector<Point<T>> &points, tvec3<T> cellCenter, T radius);
};

/* ////////////////////////////////////////////
*  //////////////   DEFINITIONS   /////////////
*/ ////////////////////////////////////////////

template <class T>
bool GeneralQuadratic<T>::fitOnPoints(vector<Point<T>> &points, tvec3<T> cellCenter, T cellSize, T radius) {
	// Calculate average point from the support points
	tvec3<T> avgPoint = MPUIUtility::positionAverage(points);

	// Set up A matrix and b vector for SVD
	Matrix<T, 10, 10> A;
	Matrix<T, 10, 1> b;

	A.setZero();
	b.setZero();

	/* /////////////////// PROCESS AUXILARY POINTS /////////////////// */

	// Setup auxilary points
	vector<tvec3<T> > auxP = MPUIUtility::generateAuxilaryPoints(cellCenter, cellSize);

	// Will remain true if no valid auxilary point is found
	bool auxEmpty = true;

	for (auto apIt = auxP.begin(); apIt != auxP.end(); apIt++) {
		// Sort ascending by distance from the auxilary point
		sort(points.begin(), points.end(), 
		[apIt](const Point<T> &pA, const Point<T> &pB) -> bool {
			tvec3<T> a = pA.position - (*apIt);
			tvec3<T> b = pB.position - (*apIt);
			return dot(a, a) < dot(b, b);
		});

		// Setup the initial sign
		bool sign = 0 >= dot(points.begin()->normal, ((*apIt) - points.begin()->position));
		bool valid = true;

		// Cumulative distance
		T dMean = 0;

		// Check auxN nearest points
		for (auto pointsIt = points.begin(); pointsIt != points.begin() + this->auxN; pointsIt++) {
			T d = dot(pointsIt->normal, ((*apIt) - pointsIt->position));

			// Check if the sign matches
			if (sign != (0 >= d) || d == 0) {
				valid = false;
				break;
			}

			dMean += d;
		}

		// Check if the auxilary point is valid. If not discard it.
		if (!valid) {
			continue;
		}

		// Found a valid auxilary point
		auxEmpty = false;

		// Calculate mean cumulative distance
		dMean /= this->auxN;

		// Auxilary point relative to the average point
		tvec3<T> avgRel = (*apIt) - avgPoint;

		// Calculate the components
		T components[10];

		components[0] = this->auxWeight;
		components[1] = this->auxWeight * avgRel.x;
		components[2] = this->auxWeight * avgRel.y;
		components[3] = this->auxWeight * avgRel.z;
		components[4] = this->auxWeight * avgRel.x * avgRel.x;
		components[5] = this->auxWeight * avgRel.y * avgRel.y;
		components[6] = this->auxWeight * avgRel.z * avgRel.z;
		components[7] = this->auxWeight * avgRel.x * avgRel.y;
		components[8] = this->auxWeight * avgRel.y * avgRel.z;
		components[9] = this->auxWeight * avgRel.z * avgRel.x;

		// Write to A matrix and B vector
		for (int i = 0; i < 10; i++) {
			for (int j = i; j < 10; j++) {
				A(i, j) += components[i] * components[j];
			}

			b(i) += components[i] * dMean * this->auxWeight;
		}
	}

	// Notify that the quadric could not be fitted due to the lack of valid auxilary points
	if (auxEmpty) {
		this->initialized = false;
		return false;
	}



	/* /////////////////// PROCESS SELECTED POINTS /////////////////// */
	T totalWeight = 0;

	// System rows for selected points
	Matrix<T, 10, 10> selPA;
	Matrix<T, 10, 1> selPb;

	selPA.setZero();
	selPb.setZero();

	for (auto pIt = points.begin(); pIt != points.end(); pIt++) {
		tvec3<T> position = pIt->position;

		// Calculate the weight for the point
		T weight = MPUIUtility::weight(position, cellCenter, radius);
		totalWeight += weight;

		// Point relative to the average point
		tvec3<T> avgRel = position - avgPoint;

		// Calculate the components
		T components[10];

		components[0] = weight;
		components[1] = weight * avgRel.x;
		components[2] = weight * avgRel.y;
		components[3] = weight * avgRel.z;
		components[4] = weight * avgRel.x * avgRel.x;
		components[5] = weight * avgRel.y * avgRel.y;
		components[6] = weight * avgRel.z * avgRel.z;
		components[7] = weight * avgRel.x * avgRel.y;
		components[8] = weight * avgRel.y * avgRel.z;
		components[9] = weight * avgRel.z * avgRel.x;

		for (int i = 0; i < 10; i++) {
			for (int j = i; j < 10; j++)
				selPA(i, j) += components[i] * components[j];
		}
	}

	// Divide the selected point matrix with total weight
	selPA /= totalWeight;

	// Combine the matrix
	A += selPA;

	// A is upper triangular matrix. We mirror it and make it symetrical
	for (int i = 1; i < 10; i++) {
		for (int j = 0; j < i; j++) {
			A(i, j) = A(j, i);
		}
	}


	/* /////////////////// SVD /////////////////// */
	
	JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);

	// Check if there is any singular value > 0
	Matrix<T, Dynamic, 1> sv = svd.singularValues();
	bool validSV = false;
	
	for (int i = 0; i < sv.rows(); i++) {
		if (sv(i, 0) > 1E-12) {
			validSV = true;
			break;
		}
	}

	if (!validSV) {
		this->initialized = false;
		return false;
	}

	// Setting the threshold for solving (if |singular value| < threshold * |max singular value| then the singular value is considered 0) 
	svd.setThreshold(this->svThresh);

	Matrix<T, Dynamic, 1> x = svd.solve(b);

	this->c[0] = x(4);	// cxx
	this->c[3] = x(5);	// cyy
	this->c[5] = x(6);	// czz

	this->c[1] = x(7);	// cxy
	this->c[4] = x(8);	// cyz
	this->c[2] = x(9);	// czx

	this->c[6] = x(1) - this->c[1] * avgPoint.y - this->c[2] * avgPoint.z - 2.0f * this->c[0] * avgPoint.x;
	this->c[7] = x(2) - this->c[4] * avgPoint.z - this->c[1] * avgPoint.x - 2.0f * this->c[3] * avgPoint.y;
	this->c[8] = x(3) - this->c[2] * avgPoint.x - this->c[4] * avgPoint.y - 2.0f * this->c[5] * avgPoint.z;

	this->c[9] = x(0) - x(1) * avgPoint.x - x(2) * avgPoint.y - x(3) * avgPoint.z
		+ this->c[1] * avgPoint.x  * avgPoint.y + this->c[4] * avgPoint.y  * avgPoint.z + this->c[2] * avgPoint.z  * avgPoint.x
		+ this->c[0] * avgPoint.x  * avgPoint.x + this->c[3] * avgPoint.y  * avgPoint.y + this->c[5] * avgPoint.z  * avgPoint.z;

	this->initialized = true;

	// Calculate the max error of the approximation
	calculateError(points, cellCenter, radius);

	return true;
}

template <class T>
void GeneralQuadratic<T>::calculateError(vector<Point<T>> &points, tvec3<T> cellCenter, T radius) {
	this->maxError = 0;

	T rSquared = radius * radius;

	for (auto itP = points.begin(); itP != points.end(); itP++) {
		tvec3<T> relPos = itP->position - cellCenter;

		if (rSquared < dot(relPos, relPos)) {
			continue;
		}

		// Caluclate value and gradient
		T fVal = abs(funValue(itP->position));
		tvec3<T> grad = gradient(itP->position);

		T l = sqrt(dot(grad, grad));
		T error;
			// Compute error for this point
			if (l != 0) {
				error = fVal / sqrt(dot(grad, grad));
			}
			else {
				error = fVal;
			}

			// Keep only the maximal error
			if (error > this->maxError) {
				this->maxError = error;
			}
	}
}


template <class T>
tvec3<T> GeneralQuadratic<T>::gradient(tvec3<T> position) {
	return tvec3<T>(
		this->c[6] + 2.0f * this->c[0] * position.x + this->c[1] * position.y + this->c[2] * position.z,
		this->c[7] + 2.0f * this->c[3] * position.y + this->c[4] * position.z + this->c[1] * position.x,
		this->c[8] + 2.0f * this->c[5] * position.z + this->c[2] * position.x + this->c[4] * position.y
		);
}


template <class T>
T GeneralQuadratic<T>::funValue(tvec3<T> position) {
	// TODO:: Check if this is working correctly.
	return  c[9] +
		(position.x * c[0] + position.y * c[1] + position.z * c[2] + c[6]) * position.x +
		(position.x * c[1] + position.y * c[3] + position.z * c[4] + c[7]) * position.y +
		(position.x * c[2] + position.y * c[4] + position.z * c[5] + c[8]) * position.z;
}

