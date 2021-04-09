#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double pi = 3.14159265358979323846;
double e = 2.7182818284590452354;

#ifdef _WIN32
// If compiling on Windows turn off warnings about fopen & fscanf
// Should be noted this is like swatting a fly with a nuke
// as it turns off a lot of Windows specific warnings
#define _CRT_SECURE_NO_DEPRECATE
#endif

// -----------------------------------------------------
// Implementing some simple complex number functionality

typedef struct Complex {
    double r; // Real part;
    double i; // Imaginary part;
} Complex;

// Calculate square of modulus for complex numbers.
double cMod(const Complex z) { return sqrt(z.r * z.r + z.i * z.i); }
double cModPtr(const Complex* z) { return cMod(*z); }

void printComplex(const Complex* z)
{
    printf("(%lf + %lfi)", z->r, z->i);
}

// End of complex number functions
// ------------------------------------------------------
// Simple struct for returning data in 3j+

typedef struct Measurement {
    size_t n;
    double time;
    Complex z;
} Measurement;

// End of Measurement struct/functions
// -------------------------------------------------------
// General functionality:
// Functions which facilitate completing the questions

// Generic write error message, then exit execution
void errorExit(const char* err_msg)
{
    fprintf(stderr, err_msg);
    exit(-1);
}

// Compares two size_t
int compareSize_t(const void* a, const void* b)
{
    return (*(size_t*)a - *(size_t*)b);
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

// Returns 1 (true) if x is in ls
// else returns 0
// Assumes ls is a sorted list
int checkIdx(const size_t* ls, const size_t sz, const size_t x)
{
    size_t i;
    for (i = 0; i < sz; ++i) {
	if (ls[i] >= x) {
	    // if ls[i]==x, x is in the list
	    if (ls[i] == x) {
		return 1;
	    } else {
		// if ls[i] > x, it isn't in the list
		return 0;
	    }
	}
    }
}

int checkIdx_binarysearch(const size_t* ls, const size_t sz, const size_t x)
{
    size_t current_idx, left_idx = 0, right_idx = sz - 1;

    // Quick fix to prevent segfault caused when sz = 1 (sz = number of elements in ls)
    if (sz == 1) {
	if (ls[0] == x)
	    return 1;
	else
	    return 0;
    }

    while (left_idx <= right_idx) {
	current_idx = floor((left_idx + right_idx) / 2);
	if (ls[current_idx] < x) {
	    left_idx = current_idx + 1;
	} else if (ls[current_idx] > x) {
	    right_idx = current_idx - 1;
	} else {
	    return 1;
	}
    }
    return 0;
}

// Produces range of N long doubles from start -> end exclusive. Simple version of numpy.linspace
double* linspaceD(const double start, const double end, const size_t N)
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
    size_t i;
    for (i = 0; i < N; i++) {
	arr[i] = val;
	val += increment;
    }

    return arr;
}

// Same as linspaceD, except values are generated using provided function, func.
Complex* linspaceComplex(Complex (*f)(double), const double* times, const size_t N)
{
    Complex* arr = (Complex*)malloc(N * sizeof(Complex));
    if (!arr) {
	errorExit("\n<linespaceComplex> malloc failed.\n");
    }
    size_t i;
    for (i = 0; i < N; ++i) {
	arr[i] = f(times[i]);
    }
    return arr;
}

// Discrete Fourier Transform taking array of samples
Complex* DFT(const Complex* samples, const size_t N)
{
    // Alloc memory for resulting array & declare vars
    Complex* arr = (Complex*)malloc(N * sizeof(Complex));
    if (!arr) {
	errorExit("\n<DFT> malloc failed.\n");
    }

    // H_n & h_k keep track of which values are currently in use
    // theta & theta_k store value of exponent h_k(t_k).exp(-2.pi.n.k/N) for use in Euler's formula
    Complex H_n, h_k;
    double theta, theta_k;
    size_t n, k;

    // For every values H_n, sum contributions of all h_k then add to resulting array
    for (n = 0; n < N; n++) {
	// Precalculate -2 * pi * n / N to prevent repetition in nested for loop
	theta = -2. * pi * n / N;

	// initialize H_n
	H_n.r = 0.;
	H_n.i = 0.;

	for (k = 0; k < N; k++) {
	    // Adjust value of theta for current sample
	    theta_k = theta * k;

	    // Retrieve k'th sample
	    h_k = samples[k];

	    // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
	    // h_k(t_k).exp(-2.pi.n.k/N):
	    H_n.r += (h_k.r * cos(theta_k) - h_k.i * sin(theta_k));
	    H_n.i += (h_k.i * cos(theta_k) + h_k.r * sin(theta_k));
	}

	arr[n] = H_n;
    }

    return arr;
}

// Placeholder function with similar function to 'pass' statement in Python
// Does nothing by design, e.g. for use in if/else or switch statements
void pass() { }

// Calculates IFT of N provided samples.
// samples --> Array of pointers to Complex values to transform.
// N --> Number of elements in samples.
// skip_n --> Array of indexes to skip or include, depending on what skip is set to
// sz --> Number of elements in skip_n
Complex* IFT(const Complex* samples, size_t N, size_t* skip_n, size_t sz)
{
    printf("\n\n\n-------\n%ld\n-------\n\n\n", sz);
    // Alloc memory for resulting array
    Complex* arr = (Complex*)malloc(N * sizeof(Complex));
    if (!arr) {
	errorExit("\n<IFT> malloc failed.\n");
    }

    // Sorting skip_n to allow checkIdx to be more efficient
    // qsort((void*)skip_n, sz, sizeof(size_t), compareSize_t);

    // H_n & h_k keep track of which values are currently in use
    // theta & theta_k store value of exponent h_k(t_k).exp(-2.pi.n.k/N) for use in Euler's formula
    Complex H_n, h_k;
    double theta, theta_k;
    size_t n, k;

    // For every values H_n, sum contributions of all h_k then add to resulting array
    for (k = 0; k < N; k++) {
	// Precalculate 2 * pi * n / N to prevent repetition in nested for loop
	theta = 2. * pi * k / N;

	// initialize h_k
	h_k.r = 0.;
	h_k.i = 0.;

	for (n = 0; n < N; n++) {
	    if (checkIdx_binarysearch(skip_n, sz, n)) {
		pass();
	    } else {
		// Adjust value of theta for current sample
		theta_k = theta * n;

		// Retrieve k'th sample
		H_n = samples[n];

		// (a + bi)(c + di) = (ac - bd) + (ad + bc)i
		// h_k(t_k).exp(-2.pi.n.k/N):
		h_k.r += (H_n.r * cos(theta_k) - H_n.i * sin(theta_k));
		h_k.i += (H_n.r * sin(theta_k) + H_n.i * cos(theta_k));
	    }
	}

	h_k.r = (double)h_k.r / N;
	h_k.i = (double)h_k.i / N;

	arr[k] = h_k;
    }

    return arr;
}

// Write Complex data to specified file
void writeComplex(const double* time_data, const Complex* z_data, const size_t N, const char* filename)
{
    // Open file, check for failure
    FILE* fp = fopen(filename, "w");
    if (!fp) {
	errorExit("\n<writeComplex> Failed to open file.\n");
    }

    // Print column titles
    fprintf(fp, "time,real,imag\n");

    // Write time & Complex data
    size_t i;
    for (i = 0; i < N; ++i) {
	fprintf(fp, "%lf,%lf,%lf\n", time_data[i], z_data[i].r, z_data[i].i);
    }

    // Close file
    fclose(fp);
}

// End of general functions
// -----------------------------------------------------
// Question answers

// Implementation for 3b
Complex** q_3b(double* times, size_t N)
{
    // Allocating memory for arrays to store results in
    Complex* h1_vals = (Complex*)malloc(N * sizeof(Complex));
    Complex* h2_vals = (Complex*)malloc(N * sizeof(Complex));
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));

    if (!h1_vals || !h2_vals || !results) {
	errorExit("\n<q_3b> malloc failed.\n");
    }

    // Generate values of h_1(t) & h_2(t)
    h1_vals = linspaceComplex(h_1, times, N);
    h2_vals = linspaceComplex(h_2, times, N);

    // Write to file
    writeComplex(times, h1_vals, N, "h1.txt");
    writeComplex(times, h2_vals, N, "h2.txt");

    // Add to results array so they can be returned and used later

    results[0] = h1_vals;
    results[1] = h2_vals;

    return results;
}

// Questions 3.d & 3.e:
// 	- 3.d: Perform DFT on h1(t) & h2(t) --> H1(w) & H2(w)
// 	       typeof(Hn(w)) = Complex*
// 	       Calulated using DFT_functions
// 	- 3.e: Print results to screen
// Also returns array of {H1(w), H2(w)} for use in part 3.f
Complex** q_3de(const double* times, Complex** samples, const size_t N)
{
    size_t i;
    // Allocating memory for results; will store {H1, H2}
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
    if (!results) {
	errorExit("\n<q_3de> malloc failed.\n");
    }

    // Calculating H1 & H2 using DFT_function
    Complex* H1 = DFT(samples[0], N);
    Complex* H2 = DFT(samples[1], N);

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

// Applies Inverse Fourier Transform (IFT) to provided sets of H_x in samples.
// In this case samples will provided by result of q_3de.
Complex** q_3f(Complex** samples, size_t N)
{
    printf("\t1\n");
    // Unpack H1 & H2 from samples.
    Complex* H1 = samples[0];
    Complex* H2 = samples[1];

    printf("\t2\n");
    // Arrays containing the indexes to skip in H1 & H2
    size_t h1_prime_skip[1] = { 1 };
    size_t h2_prime_skip[1] = { 0 };

    printf("\t3\n");
    // Calculate Inverse Fourier Transform of H1 & H2 respectively.
    Complex* h1_prime = IFT(H1, N, h1_prime_skip, 1);
    Complex* h2_prime = IFT(H2, N, h2_prime_skip, 1);

    printf("\t4\n");
    // Allocate memory for results, check for fail.
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
    if (!results) {
	errorExit("\n<q_3f> malloc failed.\n");
    }

    printf("\t5\n");
    // Pack h1_prime & h2_prime into results so it can be returned for q_3f
    results[0] = h1_prime;
    results[1] = h2_prime;

    printf("\t5\n");
    return results;
}

void q_3g(Complex** data, double* times, size_t N)
{
    // Write data to respective files
    writeComplex(times, data[0], N, "inv_1.txt");
    writeComplex(times, data[1], N, "inv_2.txt");
}

// Read data in from h3.txt, return as array of Complexes
Measurement* q_3i(const char* filename, size_t N)
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
    int n;
    double time;
    Complex z;
    size_t sz = 0;

    double real = 0.;
    double imag = 0.;

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
Measurement* q_3j(const Measurement* data, const size_t N)
{
    Complex* z_arr = (Complex*)malloc(N * sizeof(Complex));
    if (!data) {
	errorExit("\n<q_3j> malloc failed.\n");
    }

    // Copying data into array
    // DFT_samples is then applied to it
    size_t i;
    for (i = 0; i < N; ++i) {
	z_arr[i] = data[i].z;
    }

    // Applying DFT to data with DFT_samples
    Complex* DFT_data = DFT(z_arr, N);

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

// Applies IFT to 4 terms of H3 with largest amplitude
Complex* q_3k(Measurement* samples, const size_t N, size_t n_largest_vals)
{
    // Allocating memory for transformed version of data
    Complex* result = (Complex*)malloc(N * sizeof(Complex));
    if (!result) {
	errorExit("\n<q_3k> malloc failed.\n");
    }

    // Sorting a copy of the samples
    // Allocating memory for copy of array to be sorted
    Measurement* samples_sorted = (Measurement*)malloc(N * sizeof(Measurement));

    // Copying raw bits of samples --> samples_sorted
    // memcpy returns ptr --> samples_sorted, however as it already copied
    // the bits of samples over,  this is unnecessary
    memcpy((void*)samples_sorted, (void*)samples, N * sizeof(Measurement));

    // Sorting samples_sorted based on magnitude of Complex numbers stored within
    qsort(samples_sorted, N, sizeof(Measurement), compareMeasurement);

    // Variables required for
    size_t i, n, *skip_n;

    n = N - n_largest_vals;
    if (n <= 0) {
	errorExit("\n<q_3k> Performing IFT on n <= 0 values.");
    }

    skip_n = (size_t*)malloc(n * sizeof(size_t));
    if (!skip_n) {
	errorExit("\n<q_3k> malloc failed.\n");
    }

    for (i = 0; i < N; ++i) {
	skip_n[i] = samples_sorted[i].n;
    }

    printf("\t6\n");
    // Create array of Complex vals from original samples set
    Complex* z_arr = (Complex*)malloc(N * sizeof(Complex));
    if (!z_arr) {
	errorExit("\n<q_3k> malloc failed.\n");
    }
    for (i = 0; i < N; ++i) {
	z_arr[i] = samples[i].z;
    }

    printf("\t7\n");
    result = IFT(z_arr, N, skip_n, 196);

    printf("\t8\n");
    return result;
}

void q_3l(const double* times, const Complex* data, const size_t N, const char* filename)
{
    writeComplex(times, data, N, filename);
}

// End of question answers
// -----------------------------------------------------
// main()

int main(int argc, char* argv[])
{
    // Set number of samples & generate time values
    size_t i, N = 100;
    double* times = linspaceD(0, 2 * pi, N);

    printf("\nvalue of N = %ld.\n", N);

    // Complete Q3.b, then pass the values --> array of array of Complexes
    Complex** h1_and_h2 = q_3b(times, N);

    // Complete Q3.d & Q3.e, then pass the values --> array of array of Complexes
    Complex** H1_and_H2 = q_3de(times, h1_and_h2, N);

    printf("1\n");

    Complex** prime_h1_and_h2 = q_3f(H1_and_H2, N);

    printf("2\n");

    q_3g(prime_h1_and_h2, times, N);

    printf("3\n");

    // Data provided uses N=200
    N = 200;

    Measurement* h3 = q_3i("h3.txt", N);

    printf("4\n");

    Measurement* H3 = q_3j(h3, N);

    printf("5\n");

    Complex* i_h3 = q_3k(H3, N, 4);

    printf("6\n");

    times = (double*)realloc(times, N * sizeof(double));
    if (!times) {
	errorExit("\nFailed malloc for times in main.\n");
    }

    for (i = 0; i < N; ++i) {
	times[i] = h3[i].time;
    }

    printf("7\n");

    q_3l(times, i_h3, N, "inv_3.txt");

    printf("8\n");

    /*N = 128;
    size_t n = 100;

    Complex* h1 = h1_and_h2[0];
    Complex* test = (Complex*)malloc(N * sizeof(Complex));

    for (i = 0; i < N; ++i) {
	if (i < n)
	    test[i] = h1[i];
	else {
	    Complex a;
	    a.r = 0.;
	    a.i = 0.;
	    printComplex(&a);
	    test[i] = a;
	}
    }

    Complex* H1 = DFT(test, N);

    for (i = 0; i < N; ++i) {
	printf("%ld: ", i);
	printComplex(&h1[i]);
	putchar('\t');
	printComplex(&H1[i]);
	putchar('\n');
    }*/

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
