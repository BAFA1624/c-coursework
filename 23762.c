#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

// Global variables for pi & e
double pi = 3.14159265358979323846;
double e = 2.7182818284590452354;

// Generic write error message, then exit execution
void errorExit(const char* err_msg)
{
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

Complex cConjugate(Complex* z)
{
	z->i *= -1;
	return *z;
}

void printComplex(const Complex* z)
{
	printf("(%lf + %lfi)", z->r, z->i);
}

// Compares the amplitudes of two Complex numbers
int compareComplex(const void* a, const void* b)
{
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
void writeComplex(const double* time_data, const Complex* z_data, const int N, const char* filename)
{
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
int compareMeasurement(const void* a, const void* b)
{
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
int compareInt(const void* a, const void* b)
{
	return *(int*)a - *(int*)b;
}

// Implementation of h1 & h2
Complex h_1(const double time)
{
	// Calculate h1 using Euler's formula.
	Complex result = {
		.r = cos(time) + cos(5 * time),
		.i = sin(time) + sin(5 * time)
	};
	return result;
}

Complex h_2(const double time)
{
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
int checkIdx(const int* ls, const int sz, const int x)
{
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
double* linspaceD(const double start, const double end, const int N)
{
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
Complex* linspaceComplex(Complex (*f)(double), const double* samples, const int N)
{
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
Complex* DFT(const Complex* samples, const int N, const int* skip_n, const int sz, const bool IDFT)
{
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
Complex* IDFT(const Complex* samples, const int N, const int* skip_n, const int sz)
{
	// Apply DFT to samples, setting IDFT param as true and passing skip_n and sz on to DFT
	Complex* result = DFT(samples, N, skip_n, sz, true);

	// Divide all members of result by N.
	for (int i = 0; i < N; ++i) {
		result[i].r /= N;
		result[i].i /= N;
	}

	return result;
}

// End of general functions
// -----------------------------------------------------
// Question answers

// Implementation for 3b
Complex** q_3b(const double* times, const int N)
{
	// Allocating memory for arrays to store results in
	Complex** results = (Complex**)malloc(2 * sizeof(Complex*));

	if (!results) {
		errorExit("\n<q_3b> malloc failed.\n");
	}

	// Generate values of h_1(t) & h_2(t)
	Complex* h1_vals = linspaceComplex(h_1, times, N);
	Complex* h2_vals = linspaceComplex(h_2, times, N);

	// Write to file
	writeComplex(times, h1_vals, N, "h1.txt");
	writeComplex(times, h2_vals, N, "h2.txt");

	// Add to results array so they can be returned and used later
	results[0] = h1_vals;
	results[1] = h2_vals;

	return results;
}

// Questions 3.d & 3.e:
//	- 3.d: Perform DFT on h1(t) & h2(t) --> H1(w) & H2(w)
//	       typeof(Hn(w)) = Complex*
//	- 3.e: Print results to screen
// Also returns array of {H1(w), H2(w)} for use in part 3.f
Complex** q_3de(const double* times, Complex** samples, const int N)
{
	int i;
	// Allocating memory for results; will store {H1, H2}
	Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
	if (!results) {
		errorExit("\n<q_3de> malloc failed.\n");
	}

	// Calculating H1 & H2 using DFT
	Complex* H1 = DFT(samples[0], N, NULL, 0, false);
	Complex* H2 = DFT(samples[1], N, NULL, 0, false);

	// Add H1 & H2 --> results
	results[0] = H1;
	results[1] = H2;

	// Printing H1 & H2 using custom print_Complex function
	printf("\nPrinting H1:\n");
	for (i = 0; i < N; ++i) {
		printf("t = %lf\tH1 = ", times[i]);
		printComplex(&H1[i]);
		putchar('\n');
	}
	printf("\nPrinting H2:\n");
	for (i = 0; i < N; ++i) {
		printf("t = %lf\tH2 = ", times[i]);
		printComplex(&H2[i]);
		putchar('\n');
	}

	return results;
}

// Applies Inverse Fourier Transform (IDFT) to provided sets of H_x in samples.
// In this case samples will provided by result of q_3de.
Complex** q_3f(Complex** samples, int N)
{
	// Arrays containing the indexes to skip in H1 & H2
	int* h1_prime_skip = (int*)malloc(sizeof(int));
	int* h2_prime_skip = (int*)malloc(sizeof(int));

	h1_prime_skip[0] = 1;
	h2_prime_skip[0] = 0;

	// Calculate Inverse Fourier Transform of H1 & H2 respectively.
	Complex* h1_prime = IDFT(samples[0], N, h1_prime_skip, 1);
	Complex* h2_prime = IDFT(samples[1], N, h2_prime_skip, 1);

	// Freeing h1_pr... & h2_pr... to prevent memory leaks.
	free(h1_prime_skip);
	free(h2_prime_skip);

	// Allocate memory for results, check for fail.
	Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
	if (!results) {
		errorExit("\n<q_3f> malloc failed.\n");
	}

	// Pack h1_prime & h2_prime into results so it can be returned for q_3f
	results[0] = h1_prime;
	results[1] = h2_prime;

	return results;
}

void q_3g(Complex** data, double* times, int N)
{
	// Write data to respective files
	writeComplex(times, data[0], N, "inv_1.txt");
	writeComplex(times, data[1], N, "inv_2.txt");
}

// Read data in from h3.txt, return as array of Complexes
Measurement* q_3i(const char* filename, int N)
{
	FILE* fp = fopen(filename, "r");
	if (!fp) {
		errorExit("\n<q_3h> Failed to open h3.txt.\n");
	}

	// Allocating vars to keep track of array size and arrays for data
	Measurement* results = (Measurement*)malloc(N * sizeof(Measurement));
	if (!results) {
		errorExit("\n<q_3h> malloc failed.\n");
	}

	int n, sz = 0;
	double time = 0., real = 0., imag = 0.;
	Complex z;

	// Use fscanf to read from file
	while (fscanf(fp, "%d, %lf, %lf, %lf", &n, &time, &real, &imag) != EOF && sz < N) {
		// Adding values to results
		z.r = real;
		z.i = imag;

		results[sz].n = n;
		results[sz].time = time;
		results[sz].z = z;

		sz++;
	}

	fclose(fp);

	return results;
}

// Apply DFT to h3
Measurement* q_3j(const Measurement* data, const int N)
{
	Complex* z_arr = (Complex*)malloc(N * sizeof(Complex));
	if (!z_arr) {
		errorExit("\n<q_3j> malloc failed.\n");
	}

	// Copying data into array
	// DFT is then applied to it
	int i;
	for (i = 0; i < N; ++i) {
		z_arr[i] = data[i].z;
	}

	// Applying DFT to data with DFT_samples
	Complex* DFT_data = DFT(z_arr, N, NULL, 0, false);

	Measurement* results = (Measurement*)malloc(N * sizeof(Measurement));
	if (!results) {
		errorExit("<q_3j> \nFailed malloc for results in q_3j.\n");
	}

	// Copying transformed data into new array.
	for (i = 0; i < N; ++i) {
		results[i].n = data[i].n;
		results[i].time = data[i].time;
		results[i].z = DFT_data[i];
	}

	// Freeing memory not being passed out of the scope
	free(z_arr);
	free(DFT_data);

	return results;
}

// Applies IDFT to 4 terms of H3 with largest amplitude
Complex* q_3k(Measurement* samples, const int N, int n_largest_vals)
{
	// Sorting samples_sorted based on magnitude of Complex numbers stored within
	qsort(samples, N, sizeof(Measurement), compareMeasurement);

	// Variables required for
	int i, n;
	int* skip_n;

	n = N - n_largest_vals;
	if (n <= 0 || n > N) {
		printf("n = %d\n", n);
		errorExit("\n<q_3k> Performing IDFT on invalid n.");
	}

	// Allocate memory for values to skip
	skip_n = (int*)malloc(n * sizeof(int));
	if (!skip_n) {
		errorExit("\n<q_3k> malloc failed.\n");
	}

	// Copy up to n-th value into skip_n
	for (i = 0; i < n; ++i) {
		skip_n[i] = samples[i].n;
	}

	// Resort back into order by index.
	qsort((void*)samples, N, sizeof(Measurement), compareInt);

	// Create array of Complex vals from original samples set
	Complex* z_arr = (Complex*)malloc(N * sizeof(Complex));
	if (!z_arr) {
		errorExit("\n<q_3k> malloc failed.\n");
	}

	// Copy Complex data to perform IDFT on.
	for (i = 0; i < N; ++i) {
		z_arr[i] = samples[i].z;
	}

	Complex* result = IDFT(z_arr, N, skip_n, n);

	// Free memory to prevent leak.
	free(skip_n);
	free(z_arr);

	return result;
}

void q_3l(const double* times, const Complex* data, const int N, const char* filename)
{
	writeComplex(times, data, N, filename);
}

// End of question answers
// -----------------------------------------------------
// main()

int main()
{
	// Set number of samples & generate time values
	int i, N = 100;
	double* times = linspaceD(0, 2 * pi, N);

	// Complete Q3.b, then pass the values --> array of array of Complexes
	Complex** h1_and_h2 = q_3b(times, N);

	// Complete Q3.d & Q3.e, then pass the values --> array of array of Complexes
	Complex** H1_and_H2 = q_3de(times, h1_and_h2, N);

	// Complex Q3.f, then pass the values --> array of array of Complexes
	Complex** prime_h1_and_h2 = q_3f(H1_and_H2, N);

	// Complex Q3.g
	q_3g(prime_h1_and_h2, times, N);

	// Data provided uses N=200
	N = 200;

	Measurement* h3 = q_3i("h3.txt", N);

	Measurement* H3 = q_3j(h3, N);

	Complex* i_h3 = q_3k(H3, N, 4);

	times = (double*)realloc(times, N * sizeof(double));
	if (!times) {
		errorExit("\nFailed malloc for times in main.\n");
	}

	for (i = 0; i < N; ++i) {
		times[i] = h3[i].time;
	}

	q_3l(times, i_h3, N, "inv_3.txt");

	// Freeing memory
	free(times);
	free(h1_and_h2[0]);
	free(h1_and_h2[1]);
	free(h1_and_h2);
	free(H1_and_H2[0]);
	free(H1_and_H2[1]);
	free(H1_and_H2);
	free(prime_h1_and_h2[0]);
	free(prime_h1_and_h2[1]);
	free(prime_h1_and_h2);
	free(h3);
	free(H3);
	free(i_h3);

	return 0;
}
