#include "DFT.h"

int main() {
	int64_t i = 0, sz = 0, N = 100;
	double real = 0, imag = 0;

	char* filename = "/home/ben/CLionProjects/FourierTransforms/cmake-build-debug/f.txt";
	FILE* fp = fopen(filename, "r");
	if (!fp) {
		errorExit("\nFailed to open file.\n");
	}

	Complex* z_arr = (Complex*)malloc(N * sizeof(Complex));
	if (!z_arr) {
		errorExit("\nFailed malloc.\n");
	}

	for (i = 0; i < N; ++i) {
		// printComplex(&z_arr[i]);
	}

	while (fscanf(fp, "%lf,%lf", &real, &imag) != EOF && sz < N) {
		z_arr[sz].r = real;
		z_arr[sz].i = imag;
		sz++;
	}
	fclose(fp);

	Complex* z2_arr = DFT(z_arr, N, NULL, 0, false);

	for (i = 0; i < N; ++i) {
		printComplex(&z2_arr[i]);
	}

	double* tvals = linspaceD(0, 2 * pi, N);
	writeComplex(tvals, z2_arr, N, "/home/ben/CLionProjects/FourierTransforms/cmake-build-debug/F.txt");
}
