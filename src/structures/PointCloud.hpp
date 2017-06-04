#pragma once

#include "Point.hpp"
#include <vector>
#include <glm/glm.hpp>
#include "Volume.hpp"
#include "AABB.hpp"

using namespace std;
using namespace glm;

/*
* CLASS PointCloud
*/
template <class PointT>
class PointCloud {

public:
	// List of point cloud points
	vector<Point<PointT>> points;
	// Bounding box
	AABB<PointT> aabb;

public:
	PointCloud() {};
	~PointCloud() {};

	/*
	* Converts volumetric data to point cloud representation
	*/
	template <typename T>
	void setFromVolume(Volume<T> &volume, T isovalue, bool normalize = true);

	/*
	* Functions required by KD-Tree
	*/
	inline size_t kdtree_get_point_count() const { return points.size(); }

	inline PointT kdtree_distance(const PointT *p1, const size_t idx_p2, size_t /*size*/) const
	{
		const PointT d0 = p1[0] - points[idx_p2].position.x;
		const PointT d1 = p1[1] - points[idx_p2].position.y;
		const PointT d2 = p1[2] - points[idx_p2].position.z;
		return sqrt(d0*d0 + d1*d1 + d2*d2);
	}

	inline PointT kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) {
			return points[idx].position.x;
		}
		else if (dim == 1) {
			return points[idx].position.y;
		}
		else {
			return points[idx].position.z;
		}
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }

private:
	/*
	* Masks unwanted volume artifacts
	*/
	template <typename T>
	vector<bool> &generateMask(Volume<T> &volume, T isovalue);

	/*
	* Generates points based on the input parameters and returns the number of newly added points
	*/
	int generateAndAddPoints(vector<bool> &neighbourStates, int nActive, tvec3<PointT> volPos, tvec3<PointT> gradient, vector<Point<PointT>> *points);
};

/* ////////////////////////////////////////////
*  //////////////   DEFINITIONS   /////////////
*/ ////////////////////////////////////////////
template <class PointT>
template <typename T>
void PointCloud<PointT>::setFromVolume(Volume<T> &volume, T isovalue, bool normalize) {
	vector<int> dim = volume.dimensions();

	// Used to mask insignificant voxels
	//vector<bool> mask = this->generateMask(volume, isovalue);

	// Index converter
	auto at = [dim](int x, int y, int z) { return x + dim[1] * (y + dim[2] * z); };

	// Generate point cloud
	tvec3<PointT> gradient;
	int nActive = 0;
	vector<bool> neighbourStates(6, false);

	vector<Point<PointT>> points;

	for (int i = 1; i < dim[0] - 1; i++) {
		for (int j = 1; j < dim[1] - 1; j++) {
			for (int k = 1; k < dim[2] - 1; k++) {
				gradient = tvec3<PointT>(0.0f, 0.0f, 0.0f);

				// Check if the voxel is turned off and is not masked
				if (volume[i][j][k] < isovalue /*&& !mask[at(i, j, k]*/) {
					// Check which neighbours are activated
					if (volume[i][j][k - 1] > isovalue) {
						neighbourStates[0] = true;
						gradient[2]++;
						nActive++;
					}
					if (volume[i][j][k + 1] > isovalue) {
						neighbourStates[1] = true;
						gradient[2]--;
						nActive++;
					}
					if (volume[i][j - 1][k] > isovalue) {
						neighbourStates[2] = true;
						gradient[1]++;
						nActive++;
					}
					if (volume[i][j + 1][k] > isovalue) {
						neighbourStates[3] = true;
						gradient[1]--;
						nActive++;
					}
					if (volume[i - 1][j][k] > isovalue) {
						neighbourStates[4] = true;
						gradient[0]++;
						nActive++;
					}
					if (volume[i + 1][j][k] > isovalue) {
						neighbourStates[5] = true;
						gradient[0]--;
						nActive++;
					}

					gradient = glm::normalize(gradient);

					// Adds new points based on the neighbour states (modifies points vector)
					generateAndAddPoints(neighbourStates, nActive, tvec3<PointT>(i, j, k), gradient, &points);

					// Reset neighbour vector
					neighbourStates = { false, false, false, false, false, false };
					nActive = 0;
				}
			}
		}
	}

	// Normalize the points to a cube with unit-length diagonal
	if (normalize) {
		AABB<PointT> aabb(points);
		
		tvec3<PointT> center = aabb.getCenter();
		PointT maxEdgeAABB = aabb.getLongestEdge();
		
		for (size_t i = 0; i < points.size(); i++) {
			points[i].position -= center;
			points[i].position.x /= maxEdgeAABB;
			points[i].position.y /= maxEdgeAABB;
			points[i].position.z /= maxEdgeAABB;
		}

		// Fix the aabb (actually create AAAA)
		aabb.setAABB(tvec3<PointT>(-0.5, -0.5, -0.5), tvec3<PointT>(0.5, 0.5, 0.5));

		this->aabb = aabb;
	}
	else {
		this->aabb = AABB<PointT>(points);
	}

	this->points = points;
}

template <class PointT>
int PointCloud<PointT>::generateAndAddPoints(vector<bool> &neighbourStates, int nActive, tvec3<PointT> volPos, tvec3<PointT> gradient, vector<Point<PointT>> *points) {
	int count = 0;

	switch (nActive) {
		// One neighbour is active
	case 1:
		if (neighbourStates[0]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z), gradient));
		}
		else if (neighbourStates[1]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 1.0f), gradient));
		}
		else if (neighbourStates[2]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y, volPos.z + 0.5f), gradient));
		}
		else if (neighbourStates[3]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 1.0f, volPos.z + 0.5f), gradient));
		}
		else if (neighbourStates[4]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x, volPos.y + 0.5, volPos.z + 0.5f), gradient));
		}
		else {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 1.0f, volPos.y + 0.5, volPos.z + 0.5f), gradient));
		}
		count++;
		break;

		// Two neighbours are active
	case 2:
		if (neighbourStates[0] && neighbourStates[1]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z), tvec3<PointT>(0.0f, 0.0f, 1.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 1.0f), tvec3<PointT>(0.0f, 0.0f, -1.0f)));
			count += 2;
		}
		else if (neighbourStates[2] && neighbourStates[3]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y, volPos.z + 0.5f), tvec3<PointT>(0.0f, 1.0f, 0.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 1.0f, volPos.z + 0.5f), tvec3<PointT>(0.0f, -1.0f, 0.0f)));
			count += 2;
		}
		else if (neighbourStates[4] && neighbourStates[5]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x, volPos.y + 0.5f, volPos.z + 0.5f), tvec3<PointT>(1.0f, 0.0f, 0.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 1.0f, volPos.y + 0.5f, volPos.z + 0.5f), tvec3<PointT>(-1.0f, 0.0f, 0.0f)));
			count += 2;
		}
		else {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 0.5f), gradient));
			count++;
		}
		break;

		// Three neighbours are active
	case 3:
		(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 0.5f), gradient));
		count++;
		break;

		// Four neighbours are active
	case 4:
		if (!neighbourStates[0] && !neighbourStates[1]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y, volPos.z + 0.5f), tvec3<PointT>(0.0f, 1.0f, 0.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 1.0f, volPos.z + 0.5f), tvec3<PointT>(0.0f, -1.0f, 0.0f)));

			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x, volPos.y + 0.5f, volPos.z + 0.5f), tvec3<PointT>(1.0f, 0.0f, 0.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 1.0f, volPos.y + 0.5f, volPos.z + 0.5f), tvec3<PointT>(-1.0f, 0.0f, 0.0f)));
		}
		else if (!neighbourStates[2] && !neighbourStates[3]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z), tvec3<PointT>(0.0f, 0.0f, 1.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 1.0f), tvec3<PointT>(0.0f, 0.0f, -1.0f)));

			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x, volPos.y + 0.5f, volPos.z + 0.5f), tvec3<PointT>(1.0f, 0.0f, 0.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 1.0f, volPos.y + 0.5f, volPos.z + 0.5f), tvec3<PointT>(-1.0f, 0.0f, 0.0f)));
		}
		else if (!neighbourStates[4] && !neighbourStates[5]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y, volPos.z + 0.5f), tvec3<PointT>(0.0f, 1.0f, 0.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 1.0f, volPos.z + 0.5f), tvec3<PointT>(0.0f, -1.0f, 0.0f)));

			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z), tvec3<PointT>(0.0f, 0.0f, 1.0f)));
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 1.0f), tvec3<PointT>(0.0f, 0.0f, -1.0f)));
		}
		count += 2;
		break;

		// Five neighbours are active
	case 5:
		if (!neighbourStates[0]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z + 1.0f), gradient));
		}
		else if (!neighbourStates[1]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 0.5f, volPos.z), gradient));
		}
		else if (!neighbourStates[2]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y + 1.0f, volPos.z + 0.5f), gradient));
		}
		else if (!neighbourStates[3]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 0.5f, volPos.y, volPos.z + 0.5f), gradient));
		}
		else if (!neighbourStates[4]) {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x + 1.0f, volPos.y + 0.5, volPos.z + 0.5f), gradient));
		}
		else {
			(*points).push_back(Point<PointT>(tvec3<PointT>(volPos.x, volPos.y + 0.5, volPos.z + 0.5f), gradient));
		}
		count++;
		break;
	}

	return count;
}

template <class PointT>
template <typename T>
vector<bool> &PointCloud<PointT>::generateMask(Volume<T> &volume, T isovalue) {
	vector<int> dim = volume.dimensions();

	// Used to mask insignificant voxels
	vector<bool> mask(dim[0] * dim[1] * dim[2], false);

	// Index converter
	auto at = [](int x, int y, int z) { return x + dim[1] * (y + dim[2] * z); };

	// Remove 1 voxel holes in the model
	for (int i = 1; i< dim[0] - 1; i++) {
		for (int j = 1; j < dim[1] - 1; j++) {
			for (int k = 1; k < dim[2] - 1; k++) {

				if (volume[i][j][k] > isovalue &&
					volume[i][j][k - 1] < isovalue &&
					volume[i][j][k + 1] < isovalue &&
					volume[i][j - 1][k] < isovalue &&
					volume[i][j + 1][k] < isovalue &&
					volume[i - 1][j][k] < isovalue &&
					volume[i + 1][j][k] < isovalue &&

					volume[i][j + 1][k - 1] < isovalue &&
					volume[i][j - 1][k - 1] < isovalue &&
					volume[i + 1][j][k - 1] < isovalue &&
					volume[i - 1][j][k - 1] < isovalue &&

					volume[i][j + 1][k + 1] < isovalue &&
					volume[i][j - 1][k + 1] < isovalue &&
					volume[i + 1][j][k + 1] < isovalue &&
					volume[i - 1][j][k + 1] < isovalue &&

					volume[i + 1][j - 1][k] < isovalue &&
					volume[i + 1][j + 1][k] < isovalue &&
					volume[i - 1][j - 1][k] < isovalue &&
					volume[i - 1][j + 1][k] < isovalue) {

					mask[at(i, j, k)] = true;
				}
			}
		}
	}

	return mask;
}


/* ////////////////////////////////////////////
*  ////////////   KD-TREE Adaptor   ///////////
*/ ////////////////////////////////////////////
template <typename Derived>
struct PointCloudAdaptor
{
	const Derived &obj; //!< A const ref to the data set origin

	// The constructor that sets the data set source
	PointCloudAdaptor(const Derived &obj_) : obj(obj_) { }

	// CRTP helper method
	inline const Derived& derived() const { return obj; }

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return derived().points.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t /*size*/) const {
		const float d0 = p1[0] - derived().points[idx_p2].position.x;
		const float d1 = p1[1] - derived().points[idx_p2].position.y;
		const float d2 = p1[2] - derived().points[idx_p2].position.z;
		return d0*d0 + d1*d1 + d2*d2;
	}

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim == 0) {
			return derived().points[idx].position.x;
		}
		else if (dim == 1) {
			return derived().points[idx].position.y;
		}
		else {
			return derived().points[idx].position.z;
		}
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

};