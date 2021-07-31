#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Global variables for pi & e
double pi = 3.14159265358979323846;
double e = 2.7182818284590452354;

// Generic write error message, then exit execution
void errorExit(const char* err_msg) {
	fprintf(stderr, err_msg);
	exit(-1);
}

// -----------------------------------------------------
// Implementing some simple complex number functionality

typedef struct Complex {
	double r; // Real part;
	double i; // Imaginary part;
} Complex;

// Calculate square of modulus for complex numbers.
double cMod(const Complex z) { return sqrt(z.r * z.r + z.i * z.i); }

double cModPtr(const Complex* z) { return cMod(*z); }

Complex cConjugate(Complex* z) {
	z->i *= -1;
	return *z;
}

void printComplex(const Complex* z) {
	printf("(%lf + %lfi)\n", z->r, z->i);
}

// Compares the amplitudes of two Complex numbers
int compareComplex(const void* a, const void* b) {
	// Calculate amplitude of Complex vars
	const double mod_fa = cModPtr((const Complex*)a);
	const double mod_fb = cModPtr((const Complex*)b);

	// Comparison of two Complex numbers
	// If nan values detected --> error out; else --> return appropriate value for comparison
	if (isnan(mod_fa) || isnan(mod_fb)) {
		errorExit("\n<compareComplex> nan values provided.\n");
	} else if (mod_fa > mod_fb) {
		// a > b --> 1
		return 1;
	} else if (mod_fa < mod_fb) {
		// a < b --> -1
		return -1;
	}
	// Must be equal... Return 0
	return 0;
}

// Write Complex data to specified file
void writeComplex(const double* time_data, const Complex* z_data, const int N, const char* filename) {
	// Open file, check for failure
	FILE* fp = fopen(filename, "w");
	if (!fp) {
		errorExit("\n<writeComplex> Failed to open file.\n");
	}

	// Print column titles
	fprintf(fp, "time,real,imag\n");

	// Write time & Complex data
	int i;
	for (i = 0; i < N; ++i) {
		fprintf(fp, "%lf,%lf,%lf\n", time_data[i], z_data[i].r, z_data[i].i);
	}

	// Close file
	fclose(fp);
}

// End of complex number functions
// ------------------------------------------------------
// Simple struct for returning data in 3j+

typedef struct Measurement {
	int n;
	double time;
	Complex z;
} Measurement;

// Custom comparison function for qsort with measurements
int compareMeasurement(const void* a, const void* b) {
	// a & b converted back to Measurement
	const Measurement* a2 = (Measurement*)a;
	const Measurement* b2 = (Measurement*)b;

	// Return result of comparison between the amplitudes of
	// the complex numbers stored in a2 & b2.
	return compareComplex((void*)&a2->z, (void*)&b2->z);
}

// End of Measurement struct/functions
// -------------------------------------------------------
// General functionality:
// Functions which facilitate completing the questions

// Compares two int
int compareInt(const void* a, const void* b) {
	return *(int*)a - *(int*)b;
}

// Implementation of h1 & h2
Complex h_1(const double time) {
	// Calculate h1 using Euler's formula.
	Complex result = {
		.r = cos(time) + cos(5 * time),
		.i = sin(time) + sin(5 * time)
	};
	return result;
}

Complex h_2(const double time) {
	// calculate exponent, theta, of h2
	double theta = (time - pi) * (time - pi) / 2;
	// Initialize Complex number with result.
	Complex result = {
		.r = exp(theta),
		.i = 0
	};
	return result;
}

// Implementation of binary search to for a given value in a list
// ls MUST be sorted
int checkIdx(const int* ls, const int sz, const int x) {
	long current_idx, left_idx = 0, right_idx = sz - 1;

	// Prechecks:
	// If NULL passed to ls or sz = 0, empty list --> no skipped values / value not found
	if (ls == NULL || sz == 0) {
		return 0;
	}

	// Until all indexes are checked
	while (left_idx <= right_idx) {
		// Check index in the middle between left_idx & right_idx
		current_idx = floor((left_idx + right_idx) / 2.);

		if (ls[current_idx] < x) {
			// Sorted list --> if ls[current_idx] < x, then all vals before current_idx < x.
			left_idx = current_idx + 1;
		} else if (ls[current_idx] > x) {
			// Sorted list --> if ls[current_idx] > x, then all vals after current_idx > x.
			right_idx = current_idx - 1;
		} else {
			// Value found!
			return 1;
		}
	}
	// Value not found!
	return 0;
}

// Produces range of N long doubles from start -> end exclusive. Simple version of numpy.linspace
double* linspaceD(const double start, const double end, const int N) {
	// Allocate memory for result, calculate increment size for N values in given range
	double* arr = (double*)malloc(N * sizeof(double));
	if (!arr) {
		errorExit("\n<linespaceLd> malloc failed.\n");
	}

	// Calculate increment for N evenly spaced values from start -> end
	// 'val' keeps track of current value
	double increment = (double)(end - start) / N;
	double val = start;

	// Add values to arr, incrementing val every time
	int i;
	for (i = 0; i < N; i++) {
		arr[i] = val;
		val += increment;
	}

	return arr;
}

// Same as linspaceD, except values are generated using provided function, f.
Complex* linspaceComplex(Complex (*f)(double), const double* samples, const int N) {
	// Alloc memory & check for results
	Complex* arr = (Complex*)malloc(N * sizeof(Complex));
	if (!arr) {
		errorExit("\n<linespaceComplex> malloc failed.\n");
	}

	// Generate value for every sample in samples using specified function
	int i;
	for (i = 0; i < N; ++i) {
		arr[i] = f(samples[i]);
	}
	return arr;
}

// Discrete Fourier Transform taking array of samples
Complex* DFT(const Complex* samples, const int N, const int* skip_n, const int sz, const bool IDFT) {
	// Alloc memory for resulting array & declare vars
	Complex* arr = (Complex*)malloc(N * sizeof(Complex));
	if (!arr) {
		errorExit("\n<DFT> malloc failed.\n");
	}

	// Sorting skip_n to ensure binary search in checkIdx works. Only necessary if skip_n exists
	if (skip_n != NULL && sz > 1) {
		qsort((void*)skip_n, sz, sizeof(int), compareInt);
	}

	// H_n & h_k keep track of which values are currently in use
	// theta & theta_k store value of exponent h_k(t_k).exp(-2.pi.n.k/N) for use in Euler's formula
	Complex H_n, h_k;
	double theta_nk, theta_n, theta = 2. * pi / N;

	// If performing DFT (IDFT == false), multiply theta by -1
	if (!IDFT)
		theta *= -1;

	int n, k;

	// For every values H_n, sum contributions of all h_k then add to resulting array
	for (n = 0; n < N; n++) {
		// Calculate theta for this index n
		theta_n = theta * n;

		// initialize H_n
		H_n.r = 0.;
		H_n.i = 0.;

		for (k = 0; k < N; k++) {
			if (!checkIdx(skip_n, sz, k)) {
				// Adjust value of theta for current sample
				theta_nk = theta_n * k;

				// Retrieve k'th sample
				h_k = samples[k];

				// (a + bi)(c + di) = (ac - bd) + (ad + bc)i
				// h_k(t_k).exp(-2.pi.n.k/N):
				H_n.r += (h_k.r * cos(theta_nk) - h_k.i * sin(theta_nk));
				H_n.i += (h_k.r * sin(theta_nk) + h_k.i * cos(theta_nk));
			}
		}
		// Add to result array
		arr[n] = H_n;
	}

	return arr;
}

// Calculates IDFT of N provided samples.
// samples --> Array of pointers to Complex values to transform.
// N --> Number of elements in samples.
// skip_n --> Array of indexes to skip or include, depending on what skip is set to
// sz --> Number of elements in skip_n
Complex* IDFT(const Complex* samples, const int N, const int* skip_n, const int sz) {
	// Apply DFT to samples, setting IDFT param as true and passing skip_n and sz on to DFT
	Complex* result = DFT(samples, N, skip_n, sz, true);

	// Divide all members of result by N.
	for (int i = 0; i < N; ++i) {
		result[i].r /= N;
		result[i].i /= N;
	}

	return result;
}
