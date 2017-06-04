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
class BivariateQuadratic : public Implicit<T> {
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
	BivariateQuadratic(size_t auxN = 6, T svThresh = 1E-5) : initialized(false), auxN(auxN), svThresh(svThresh) {};
	~BivariateQuadratic() {};

	bool fitOnPoints(vector<Point<T>> &points, tvec3<T> cellCenter, T radius, T errRadius, tvec3<T> avgPoint, tvec3<T> avgNormal);

	inline tvec3<T> gradient(tvec3<T> position);

	inline T funValue(tvec3<T> position);

	inline T getApproximationError() { return this->maxError; }

private:
	void calculateError(vector<Point<T>> &points, tvec3<T> cellCenter, T radius);
};

/* ////////////////////////////////////////////
*  //////////////   DEFINITIONS   /////////////
*/ ////////////////////////////////////////////

template <class T>
bool BivariateQuadratic<T>::fitOnPoints(vector<Point<T> > &points, tvec3<T> cellCenter, T radius, T errRadius, tvec3<T> avgPoint, tvec3<T> avgNormal) {

	// Calculate base for vector space U (u, v, w) where positive w is pointing towards normal
	tvec3<T> transA;
	T l;

	if (abs(avgNormal.x) < abs(avgNormal.y)) {
		l = sqrt(avgNormal.y * avgNormal.y + avgNormal.z * avgNormal.z);
		transA = tvec3<T>(0, -(avgNormal.z / l), (avgNormal.y / l));
	}
	else {
		l = sqrt(avgNormal.x * avgNormal.x + avgNormal.z * avgNormal.z);
		transA = tvec3<T>((avgNormal.z / l), 0, -(avgNormal.x / l));
	}

	tvec3<T> transB(
		avgNormal.y * transA.z - avgNormal.z * transA.y,
		avgNormal.z * transA.x - avgNormal.x * transA.z,
		avgNormal.x * transA.y - avgNormal.y * transA.x);


	// Set up A matrix and b vector for SVD
	Matrix<T, 6, 6> A;
	Matrix<T, 6, 1> b;

	A.setZero();
	b.setZero();

	/* /////////////////// PROCESS SELECTED POINTS /////////////////// */

	for (auto pIt = points.begin(); pIt != points.end(); pIt++) {
		tvec3<T> position = pIt->position;

		// Calculate the weight for the point
		T weight = MPUIUtility::weight(position, cellCenter, radius);

		// Point relative to the average point
		tvec3<T> avgRel = position - avgPoint;

		// Transfrom the point to vector space U
		tvec3<T> posU(
			dot(transA, avgRel),
			dot(transB, avgRel),
			weight * dot(avgNormal, avgRel)
		);

		// Calculate the components
		T components[6];

		components[0] = weight;
		components[1] = weight * posU.x;
		components[2] = weight * posU.y;
		components[3] = weight * posU.x * posU.x;
		components[4] = weight * posU.x * posU.y;
		components[5] = weight * posU.y * posU.y;

		// Set A matrix coefficients and b vector
		for (int i = 0; i < 6; i++) {
			for (int j = i; j < 6; j++) {
				A(i, j) += components[i] * components[j];
			}

			b(i) += components[i] * posU.z;
		}
	}

	// A is upper triangular matrix. We mirror it and make it symetrical
	for (int i = 1; i < 6; i++) {
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

	// Solve the system
	Matrix<T, Dynamic, 1> x = svd.solve(b);

	// Set up the coefficients
	this->c[0] = -(x(3) * transA.x * transA.x + x(5) * transB.x * transB.x + x(4) * transA.x * transB.x); // cxx
	this->c[3] = -(x(3) * transA.y * transA.y + x(5) * transB.y * transB.y + x(4) * transA.y * transB.y); // cyy
	this->c[5] = -(x(3) * transA.z * transA.z + x(5) * transB.z * transB.z + x(4) * transA.z * transB.z); // czz

	this->c[1] = -2 * (x(3) * transA.x * transA.y + x(5) * transB.x * transB.y + x(4) * (transA.x * transB.y + transA.y * transB.x));	// cxy
	this->c[4] = -2 * (x(3) * transA.y * transA.z + x(5) * transB.y * transB.z + x(4) * (transA.y * transB.z + transA.z * transB.y));	// cyz
	this->c[2] = -2 * (x(3) * transA.z * transA.x + x(5) * transB.z * transB.x + x(4) * (transA.z * transB.x + transA.x * transB.z));	// czx

	T cX = avgNormal.x - x(1) * transA.x - x(2) * transB.x;
	T cY = avgNormal.y - x(1) * transA.y - x(2) * transB.y;
	T cZ = avgNormal.z - x(1) * transA.z - x(2) * transB.z;

	this->c[6] = cX - this->c[1] * avgPoint.y - this->c[2] * avgPoint.z - 2 * this->c[0] * avgPoint.x;
	this->c[7] = cY - this->c[4] * avgPoint.z - this->c[1] * avgPoint.x - 2 * this->c[3] * avgPoint.y;
	this->c[8] = cZ - this->c[2] * avgPoint.x - this->c[4] * avgPoint.y - 2 * this->c[5] * avgPoint.z;


	this->c[9] = -x(0) - cX * avgPoint.x - cY * avgPoint.y - cZ * avgPoint.z
		+ this->c[1] * avgPoint.x * avgPoint.y + this->c[4] * avgPoint.y * avgPoint.z + this->c[2] * avgPoint.z * avgPoint.x
		+ this->c[0] * avgPoint.x * avgPoint.x + this->c[3] * avgPoint.y * avgPoint.y + this->c[5] * avgPoint.z * avgPoint.z;

	this->initialized = true;

	// Calculate the max error of the approximation
	calculateError(points, cellCenter, errRadius);

	return true;
}

template <class T>
void BivariateQuadratic<T>::calculateError(vector<Point<T>> &points, tvec3<T> cellCenter, T radius) {
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
tvec3<T> BivariateQuadratic<T>::gradient(tvec3<T> position) {
	return tvec3<T>(
		this->c[6] + 2.0f * this->c[0] * position.x + this->c[1] * position.y + this->c[2] * position.z,
		this->c[7] + 2.0f * this->c[3] * position.y + this->c[4] * position.z + this->c[1] * position.x,
		this->c[8] + 2.0f * this->c[5] * position.z + this->c[2] * position.x + this->c[4] * position.y
		);
}

// cxx - c[0]
// cyy - c[3]
// czz - c[5]

// cxy - c[1]
// cyz - c[4]
// czx - c[2]

// cx - c[6]
// cy - c[7]
// cz - c[8]

// c0 - c[9]

template <class T>
T BivariateQuadratic<T>::funValue(tvec3<T> position) {
	// TODO:: Check if this is working correctly.
	return  c[9] +
		(position.x * c[0] + position.y * c[1] + position.z * c[2] + c[6]) * position.x +
		(position.x * c[1] + position.y * c[3] + position.z * c[4] + c[7]) * position.y +
		(position.x * c[2] + position.y * c[4] + position.z * c[5] + c[8]) * position.z;
	/*
	return  c[9] +
		(position.x * c[0] + position.y * c[1] + c[6]) * position.x +
		(position.y * c[3] + position.z * c[4] + c[7]) * position.y +
		(position.x * c[2] + position.z * c[5] + c[8]) * position.z;*/
}

