
#include "Point.hpp"
#include <glm/glm.hpp>

template<class T>
class AABB {
public:
	AABB();
	~AABB() {};

	// Sets AABB from the given points
	AABB(const vector<Point<T> > &points);
	
	// Extends AABB with the given points
	void extend(const vector<Point<T> > &points);

	// Sets the AABB based on the input vectors
	void setAABB(tvec3<T> min, tvec3<T> max) { this->minP = min; this->maxP = max; }

	// Fetches AABB diagonal
	glm::tvec3<T> AABB<T>::getDiagonal() const;

	// Fetches magnitude longest edge of the AABB
	T getLongestEdge() const;

	tvec3<T> getCenter() const;

	// Resets the AABB
	void setNull() { this->minP = tvec3<T>(1.0); this->maxP = tvec3<T>(-1.0); }
	bool isNull() const { return this->minP.x > this->maxP.x || this->minP.y > this->maxP.y || this->minP.z > this->maxP.z; }

	// Getters
	glm::tvec3<T> getMin() { return minP; }
	glm::tvec3<T> getMax() { return maxP; }

private:
	// AABB bounds
	tvec3<T> minP;
	tvec3<T> maxP;

	T compMax(const tvec3<T> vec) const;
};

template <class T>
AABB<T>::AABB() {
	this->setNull();
}

template <class T>
AABB<T>::AABB(const vector<Point<T> > &points) {
	this->setNull();
	this->extend(points);
}



template <class T>
void AABB<T>::extend(const vector<Point<T> > &points) {
	if (points.size() <= 0) {
		return;
	}

	auto pointsIt = points.begin();
	
	// Set initial min max
	if (this->isNull()) {
		this->minP = pointsIt->position;
		this->maxP = pointsIt->position;
		pointsIt++;
	}

	for (auto it = pointsIt; it != points.end(); it++) {
		this->minP = glm::min(this->minP, it->position);
		this->maxP = glm::max(this->maxP, it->position);
	}
}

template <class T>
glm::tvec3<T> AABB<T>::getDiagonal() const {
	if (!this->isNull()) {
		return this->maxP - this->minP;
	}
	else {
		return glm::tvec3<T>(0);
	}
}

template <class T>
glm::tvec3<T> AABB<T>::getCenter() const {
	if (!this->isNull()) {
		glm::tvec3<T> d = this->getDiagonal();
		return minP + (d * 0.5);
	}
	else {
		return glm::tvec3<T>(0.0);
	}
}



template <class T>
T AABB<T>::getLongestEdge() const {
	return this->compMax(getDiagonal());
}

template <class T>
T AABB<T>::compMax(const tvec3<T> vec) const {
	T rez = vec[0];
	for (length_t i = 1, n = vec.length(); i < n; ++i)
		rez = (rez < vec[i]) ? vec[i] : rez;
	return rez;
}