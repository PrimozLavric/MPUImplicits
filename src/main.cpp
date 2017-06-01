#include "structures\volume.hpp"
#include "structures\PointCloud.hpp"
#include  <iostream>
#include <fstream>
#include <nanoflann.hpp>
#include "structures\ImplicitOctree.hpp"

using namespace std;
using namespace nanoflann;

int main() {
	int size = 100;

	Volume<double> data(size, size, size);

	double axisMin = -10;
	double axisRange = 20;

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < size; k++) {
				// actual values
				double x = axisMin + axisRange * i / (size - 1);
				double y = axisMin + axisRange * j / (size - 1);
				double z = axisMin + axisRange * k / (size - 1);
				double value = x*x + y*y - z*z - 25;
				data[i][j][k] = value;
			}
		}
	}
	
	PointCloud<double> pc;
	pc.setFromVolume(data, 0.0, true);


	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud<double> >, PointCloud<double>, 3> kdtree;

	kdtree tree(3, pc);
	tree.buildIndex();

	double query_pt[3] = { 0.5f, 0.5f, 0.5f };
	double radius = 0.1f;

	std::vector<std::pair<size_t, double>> ret_matches;

	nanoflann::SearchParams params;
	//params.sorted = false;

	const size_t nMatches = tree.radiusSearch(&query_pt[0], radius, ret_matches, params);

	return 0;

}

