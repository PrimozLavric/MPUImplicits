#include "structures\volume.hpp"
#include "structures\PointCloud.hpp"
#include <iostream>
#include <fstream>
#include "structures\ImplicitOctree.hpp"

using namespace std;
using namespace nanoflann;

/*void test() {
	vector<Point<float> > points;

	float pointsM[20][3];
	float normalsM[20][3];
	int indicesM[20];
	float centerM[3] = { 0, 0, 0 };
	float avgM[3] = { 2.38418583e-08, 0, 0 };


	tvec3<float> center(0, 0, 0);
	float size = 1;
	float radius = 1;
	for (int i = 0; i <= 19; i++) {
		points.push_back(Point<float>(tvec3<float>((i / 19.0f) - 0.5, (i * i) / (19 * 19) - 0.5, (i * i * i) / (19 * 19 * 19) - 0.5), tvec3<float>(0, 0, 1)));
		pointsM[i][0] = (i / 19.0f) - 0.5;
		pointsM[i][1] = (i * i) / (19 * 19) - 0.5;
		pointsM[i][2] = (i * i * i) / (19 * 19 * 19) - 0.5;

		normalsM[i][0] = 0;
		normalsM[i][1] = 0;
		normalsM[i][2] = 1;

		indicesM[i] = i;
	}

	//GeneralQuadratic<float> quadric;

	//quadric.fitOnPoints(points, center, size, radius);

	BivariateQuadratic<float> bivQuadric;

	bivQuadric.fitOnPoints(points, center, radius);


	BivQuadratic(pointsM, normalsM, radius, centerM, indicesM, 20, &MPUIUtility::positionAverage(points)[0], &MPUIUtility::normalAverage(points)[0]);
}*/

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


	return 0;

}

