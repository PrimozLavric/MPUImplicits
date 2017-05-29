#include "structures\volume.h"

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

}

