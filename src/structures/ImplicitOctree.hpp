#pragma once

#include <vector>
#include <glm/glm.hpp>
#include "PointCloud.hpp"
#include <nanoflann.hpp>

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
	ImplicitOctree<T> treeRef;

	// List of cell children
	ImplicitCell<T> children[8];

	// True if the cell is a leaf
	bool isLeaf;

	// Position (min cube vector) of the octree cell
	tvec3<T> position;

	// Length of the cell side
	T size;

public:
	ImplicitCell(ImplicitOctree<T>& tree, tvec3<T> position, T size);
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
ImplicitCell<T>::ImplicitCell(ImplicitOctree<T>& tree, tvec3<T> position, T size) {
	this->treeRef = tree;
	this->position = position;
	this->size = size;
}

template <class T>
void ImplicitCell<T>::build() {
	// TODO
}

template <class T>
void ImplicitCell<T>::subdivide() {
	T halfSize = this->size / 2;

	this->children[0] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y,				position.z),			halfSize);
	this->children[1] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y,				position.z),			halfSize);
	this->children[2] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y + halfSize,	position.z),			halfSize);
	this->children[3] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y + halfSize,	position.z),			halfSize);
	this->children[4] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y,				position.z + halfSize), halfSize);
	this->children[5] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y,				position.z + halfSize), halfSize);
	this->children[6] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x,				position.y + halfSize,	position.z + halfSize), halfSize);
	this->children[7] = ImplicitCell<T>(this->treeRef, tvec3<T>(position.x + halfSize,	position.y + halfSize,	position.z + halfSize), halfSize);
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
		double supportBallR;
		double lambda; // Radius iteration step
		double minPoints;
		int maxLevel;
		T maxError;

		Parameters(double supportBallR, double lambda, double minPoints, int maxLevel, T maxError) : supportBallR(supportBallR), lambda(lambda), minPoints(minPoints), maxLevel(maxLevel), maxError(maxError) {};
	};

	// Pointcloud containing the data
	PointCloud<T> pointCloud;
	// Kdtree
	kdtree pcKdtree;
	// Root cell of the tree
	ImplicitCell root;
	// Parameters
	Parameters params;

public:
	ImplicitOctree(PointCloud<T> &pointCloud, Parameters params = Parameters(0.75, 0.1, 15, 0.001));
	~ImplicitOctree() {};

	void build();
};

//////////////////////////////////////////////////////
////////////  DEFINITIONS ImplicitOctree  ////////////
//////////////////////////////////////////////////////

template <class T>
ImplicitOctree<T>::ImplicitOctree(PointCloud<T> &pointCloud, Parameters params) {
	this->pointCloud = pointCloud;
	this->params = params;

	// Build the kdtree
	this->pcKdtree = kdtree(3, this->pointCloud);
	this->pcKdtree.buildIndex();

	// Initialize root cell
	this->root = ImplicitCell<T>(this, pointCloud.aabb.getMin(), pointCloud.aabb.getLongestEdge());
}

template <class T>
void ImplicitOctree<T>::build() {
	this->root.build();
}