#include "structures\volume.hpp"
#include "structures\PointCloud.hpp"
#include <iostream>
#include <fstream>
#include "structures\ImplicitOctree.hpp"
#include "polygonization\Bloomenthal.hpp"

using namespace std;
using namespace nanoflann;

void globalDistanceFunctionTest() {
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

	PointCloud<float> pc;
	pc.setFromVolume(data, 0.0, true);

	ImplicitOctree<float> implicitOctree(pc);
	implicitOctree.build();



	cout << "X: " << endl;
	for (int i = 0; i <= 20; i++) {
		float distance = implicitOctree.globalFunctionValue(tvec3<float>(i / 20.0 - 0.5, 0, 0));
		printf("%.6f\n", distance);
	}

	cout << "Y: " << endl;
	for (int i = 0; i <= 20; i++) {
		float distance = implicitOctree.globalFunctionValue(tvec3<float>(0, i / 20.0 - 0.5, 0));
		printf("%.6f\n", distance);
	}

	cout << "Z: " << endl;
	for (int i = 0; i <= 20; i++) {
		float distance = implicitOctree.globalFunctionValue(tvec3<float>(0, 0, i / 20.0 - 0.5));
		printf("%.6f\n", distance);
	}
}

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

	PointCloud<float> pc;
	pc.setFromVolume(data, 0.0, true);

	ImplicitOctree<float> implicitOctree(pc);
	implicitOctree.build();

	auto distLambda = [&implicitOctree](tvec3<float> x) { return implicitOctree.globalFunctionValue(x); };

	vector<Point<float> > triangulation = Bloomenthal::polygonize<float>(10, distLambda);

	ofstream myfile;
	myfile.open("out.obj");
	for (int i = 0; i < triangulation.size(); i++) {
		myfile << "v " << triangulation[i].position.x << " " << triangulation[i].position.y << " " << triangulation[i].position.z << endl;
		myfile << "n " << triangulation[i].normal.x << " " << triangulation[i].normal.y << " " << triangulation[i].normal.z << endl;
	}

	for (int i = 0; i < triangulation.size(); i += 3) {
		myfile << "f " << i+1 << " " << i+2 << " " << i+3 << endl;
	}

	
	myfile.close();

	return 0;

}

