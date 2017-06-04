#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "PointCloud.hpp"
#include <nanoflann.hpp>
#include "../implicits/Implicit.hpp"
#include "../implicits/GeneralQuadratic.hpp"
#include "../implicits/BivariateQuadratic.hpp"
#include "../implicits/UtilityFunctions.hpp"

using namespace std;
using namespace glm;
using namespace nanoflann;

// Forward declaration for usage in Implicit cell
template <class T>
class ImplicitOctree;

//////////////////////////////////////////////////////
/////////////  DECLARATION ImplicitCell  /////////////
//////////////////////////////////////////////////////

template <class T>
class ImplicitCell {
private:
	// Reference to the tree
	ImplicitOctree<T> *treeRef;

	// List of cell children
	vector<ImplicitCell<T> > children;

	// True if the cell is a leaf
	bool isLeaf;

	// Position (min cube vector) of the octree cell
	tvec3<T> position;

	// Length of the cell side
	T size;

	// Support radius
	T radius;

	// Depth of this cell
	int depth;

	// Implicit approximating points in this cell. (Only leafs should have this)
	Implicit<T> *implicit;

public:
	ImplicitCell() {};
	ImplicitCell(ImplicitOctree<T> *tree, tvec3<T> position, T size, int depth);
	~ImplicitCell() {};

	tvec3<T> getCenter() { return this->position + (this->size / 2); }

	// Subdivides the cell (creating 8 children)
	void subdivide();

	void build();
};

//////////////////////////////////////////////////////
/////////////  DEFINITIONS ImplicitCell  /////////////
//////////////////////////////////////////////////////

template <class T>
ImplicitCell<T>::ImplicitCell(ImplicitOctree<T> *tree, tvec3<T> position, T size, int depth) {
	this->treeRef = tree;
	this->position = position;
	this->size = size;
	this->depth = depth;
	this->isLeaf = true;

	this->radius = tree->params.alpha * size * sqrt(3);
}

template <class T>
void ImplicitCell<T>::build() {

	// Reference t o pointCloud points
	vector<Point<T> > pcPoints = this->treeRef->pointCloud.points;
	
	tvec3<T> cellCenter = this->getCenter();
	T query_pt[3] = { cellCenter.x, cellCenter.y, cellCenter.z };

	// List of support points
	int minPoints = this->treeRef->params.minPoints;
	vector<Point<T> > supportPoints;

	T origRadius = this->radius;

	// Fetch the points from support radius using kdtree
	// Will hold pair for each point where the first element is point index and second is the distance from that point
	vector<pair<size_t, T> > ret_matches;
	size_t nMatches = this->treeRef->pcKdtree.radiusSearch(&query_pt[0], this->radius, ret_matches, SearchParams());
	
	// Check if enough support points were found
	if (nMatches < minPoints) {
		// To few points were inside the support radius. Select minPoints nearest neighbours and set the radius accordingly
		vector<size_t> indices(minPoints);
		vector<T> distances(minPoints);

		this->treeRef->pcKdtree.knnSearch(&query_pt[0], minPoints, &indices[0], &distances[0]);

		// Set the radius to maximum distance
		this->radius = sqrt(distances[minPoints - 1]);

		// Add points to support vector points
		for (auto idxIt = indices.begin(); idxIt != indices.end(); idxIt++) {
			supportPoints.push_back(pcPoints[*idxIt]);
		}
	}
	else {
		for (auto matchIt = ret_matches.begin(); matchIt != ret_matches.end(); matchIt++) {
			supportPoints.push_back(pcPoints[matchIt->first]);
		}
	}

	// Calculate average tangent plane
	pair<tvec3<T>, tvec3<T> > tangentPlane = MPUIUtility::findAveragedTangentPlane(supportPoints, cellCenter, this->radius);

	// Check if the points should be approximated with either general quaratic or bivariate quadratic function
	bool isGeneral = false;

	for (auto itP = supportPoints.begin(); itP != supportPoints.end(); itP++) {
		if (dot(tangentPlane.second, itP->normal) < 0) {
			isGeneral = true;
			break;
		}
	}

	bool fittingSuccessful;
	if (isGeneral) {
		GeneralQuadratic<T> gq = GeneralQuadratic<T>();
		fittingSuccessful = gq.fitOnPoints(supportPoints, cellCenter, this->size, this->radius, tangentPlane.first);
		this->implicit = &gq;
		//cout << "General: " << this->implicit->getApproximationError() << "\tStatus:" << fittingSuccessful << endl;
	}
	else {
		BivariateQuadratic<T> bq = BivariateQuadratic<T>();
		fittingSuccessful = bq.fitOnPoints(supportPoints, cellCenter, this->radius, origRadius, tangentPlane.first, tangentPlane.second);
		this->implicit = &bq;
		//cout << "Bivariate: " << this->implicit->getApproximationError() << "\tStatus:" << fittingSuccessful << endl;
	}
	

	// Maximal allowed error and max depth
	T maxError = this->treeRef->params.maxError;
	int maxDepth = this->treeRef->params.maxDepth;

	//cout << this->implicit->getApproximationError() << endl;
	// If implicit fitting was unsuccessful or if the error is above the threshold subdivide this cell
	if ((!fittingSuccessful || this->implicit->getApproximationError() > maxError) && this->depth < maxDepth) {
		
		// Current approximation is not good enough.. subdivide
		this->subdivide();

		// Build sub cells
		for (int i = 0; i < 8; i++) {
			this->children[i].build();
		}
	}
	else {
		if (this->depth > this->treeRef->maxDepth) {
			this->treeRef->maxDepth = this->depth;
		}

		if (!fittingSuccessful && this->depth >= maxDepth) {
			cerr << "Maximal depth was reached and the fitting was unsuccessful!";
		}
	}
}

template <class T>
void ImplicitCell<T>::subdivide() {
	T halfSize = this->size / 2;

	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y,				position.z),			halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y,				position.z),			halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y + halfSize,	position.z),			halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y + halfSize,	position.z),			halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y,				position.z + halfSize), halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y,				position.z + halfSize), halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y + halfSize,	position.z + halfSize), halfSize, this->depth + 1));
	this->children.push_back(ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y + halfSize,	position.z + halfSize), halfSize, this->depth + 1));

	this->isLeaf = false;
}

//////////////////////////////////////////////////////
////////////  DECLARATION ImplicitOctree  ////////////
//////////////////////////////////////////////////////

template <class T>
class ImplicitOctree {

private:
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<T, PointCloud<T> >, PointCloud<T>, 3> kdtree;

public:
	struct Parameters {
		T alpha; // radius = alpha * cellDiag
		T lambda; // Radius iteration step
		int minPoints;
		int maxDepth;
		T maxError;

		Parameters(T alpha, T lambda, int minPoints, int maxDepth, T maxError) : alpha(alpha), lambda(lambda), minPoints(minPoints), maxDepth(maxDepth), maxError(maxError) {};
	};

	// Pointcloud containing the data
	PointCloud<T> pointCloud;
	// Kdtree
	kdtree pcKdtree;
	// Root cell of the tree
	ImplicitCell<T> root;
	// Parameters
	Parameters params;

	int maxDepth = 0;

public:
	ImplicitOctree(PointCloud<T> &pointCloud, Parameters params = Parameters((T) 0.75, (T) 0.1, 15, 10, (T) 0.005));
	~ImplicitOctree() {};

	void build();
};

//////////////////////////////////////////////////////
////////////  DEFINITIONS ImplicitOctree  ////////////
//////////////////////////////////////////////////////

template <class T>
ImplicitOctree<T>::ImplicitOctree(PointCloud<T> &pointCloud, Parameters params) 
	: pointCloud(pointCloud), params(params), pcKdtree(3, pointCloud) {
	root = ImplicitCell<T>(this, pointCloud.aabb.getMin(), pointCloud.aabb.getLongestEdge(), 0);
}

template <class T>
void ImplicitOctree<T>::build() {
	// Build the kdtree
	this->pcKdtree.buildIndex();

	// Start building sublevels
	this->root.build();
}