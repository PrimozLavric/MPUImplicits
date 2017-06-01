#pragma once

#include <vector>

using namespace std;

template <class T>
class Volume {

public:
	Volume(int x, int y, int z);
	~Volume() {};

	// Resizes the volume
	void resize(int x, int y, int z);

	// Initializes the volume structure with the given data
	void init(vector<T> const &data);

	// Element fetch operator
	vector<vector<T>>& operator[] (int x) {
		return this->data[x];
	}

	vector<int> const dimensions() { return this->dim; }
private:
	vector<int> dim = { -1, -1, -1 };

	// Holds the volume data
	vector<vector<vector<T>>> data;
};

template <class T>
Volume<T>::Volume(int x, int y, int z) {
	// Allocate space for the volume data
	this->data.resize(x, vector<vector<T> >(y, vector<T>(z)));

	// Store dimensions
	this->dim = { x, y, z };
}

template <class T>
void Volume<T>::resize(int x, int y, int z) {
	this->data.resize(x, vector<vector<T> >(y, vector<T>(z)));
}

template <class T>
void Volume<T>::init(vector<T> const &values) {
	auto iter = values.begin();

	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
			for (int k = 0; k < z; k++) {
				
				values[i][j][k] = *it;
				
				// Check if we have reached the end
				if (iter++ == values.end()) {
					return;
				}
			}
		}
	}
}