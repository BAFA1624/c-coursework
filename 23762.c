#include <math.h>
#include <stdio.h>
#include <stdlib.h>

long double pi = 3.14159265358979323846;
long double pi_2 = 1.57079632679489661923;
long double pi_4 = 0.78539816339744830962;
long double e = 2.7182818284590452354;

// Implementing some simple complex number functionality

typedef struct Complex {
    long double r; // Real part;
    long double i; // Imaginary part;
} Complex;

// Calculate square of modulus for complex numbers.
long double c_mod(Complex z) { return sqrtl(z.r * z.r + z.i * z.i); }
long double c_mod_ptr(Complex* z) { return c_mod(*z); }

// Calculate argument of complex number
long double c_arg(Complex z) { return atan2l(z.i, z.r); }
long double c_arg_ptr(Complex* z) { return atan2l(z->i, z->r); }

// Implementation of h1 & h2
Complex h_1(long double time)
{
    // Calculate h1 using Euler's formula.
    Complex result = {
	.r = cosl(time) + cosl(5 * time),
	.i = sinl(time) + sinl(5 * time)
    };
    return result;
}
Complex h_2(long double time)
{
    // calculate exponent, x, of h2
    long double x = (time - pi) * (time - pi) / 2;
    // Initialize Complex number with result.
    Complex result = {
	.r = powl(e, x),
	.i = 0
    };
    return result;
}

// Pops position idx from array arr of size N
long double* arr_pop_ld(long double* arr, size_t N, size_t idx)
{
    if (idx >= N) {
	fprintf(stderr, "OutOfBoundsError, attempted to pop element %ld which is out of bounds of array size %ld.\n", idx, N);
	exit(-1);
    }

    size_t i;
    for (i = idx; i < N - 1; ++i) {
	arr[i] = arr[i + 1];
    }
    arr = (long double*)realloc(arr, (N - 1) * sizeof(long double));
    if (!arr) {
	fprintf(stderr, "realloc failed in arr_pop\n");
	exit(-1);
    }
    return arr;
}

Complex* arr_pop_Complex(Complex* arr, size_t N, size_t idx)
{
    if (idx >= N) {
	fprintf(stderr, "OutOfBoundsError, attempted to pop element %ld which is out of bounds of array size %ld.\n", idx, N);
	exit(-1);
    }

    size_t i;
    for (i = idx; i < N - 1; ++i) {
	arr[i] = arr[i + 1];
    }
    arr = (Complex*)realloc(arr, (N - 1) * sizeof(Complex));
    if (!arr) {
	fprintf(stderr, "realloc failed in arr_pop\n");
	exit(-1);
    }
    return arr;
}

long double* linspace_ld(long double start, long double end, int N)
{
    long double* arr = (long double*)malloc(N * sizeof(long double));
    long double increment = (double)(end - start) / N;
    int i;
    for (i = 0; i < N; i++) {
	arr[i] = start;
	start += increment;
    }
    return arr;
}

// Discrete Fourier Transform taking sample function
Complex* DFT_functions(Complex (*func)(long double), long double* times, int N)
{
    // Takes:
    // 	- func: Ptr to function to transform.
    // 	- times: Long double array of sample times.
    // 	- N: Size of times, incidentally also number of samples.

    // Array for complex results
    Complex* results = (Complex*)malloc(N * sizeof(Complex));

    // Variables for increment and calculating Euler's formula.
    Complex H_n, h_k;
    long double theta;
    int n, k;

    for (n = 0; n < N; ++n) {
	// Initialize current value for H_n.
	H_n.r = 0.;
	H_n.i = 0.;

	for (k = 0; k < N; ++k) {
	    // Calculate values
	    h_k = func(times[k]);
	    theta = 2. * pi * n * k / N;

	    // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
	    // h_k(t_k).exp(-2.pi.n.k/N):
	    H_n.r += (h_k.r * cosl(theta) + h_k.i * sinl(theta));
	    H_n.i += (h_k.i * cosl(theta) - h_k.r * sinl(theta));
	}
	// Add to results array;
	results[n] = H_n;
    }
    return results;
}

// Discrete Fourier Transform taking array of samples
Complex* DFT_samples(Complex* samples, long double* times, int N)
{
    Complex* results = (Complex*)malloc(N * sizeof(Complex));

    Complex H_n, h_k;
    long double theta;
    int n, k;

    for (n = 0; n < N; n++) {
	H_n.r = 0.;
	H_n.i = 0.;

	for (k = 0; k < N; k++) {
	    // Calculate values
	    theta = 2. * pi * n * k / N;
	    h_k = samples[k];

	    // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
	    // h_k(t_k).exp(-2.pi.n.k/N):
	    H_n.r += (h_k.r * cosl(theta) + h_k.i * sinl(theta));
	    H_n.i += (h_k.i * cosl(theta) - h_k.r * sinl(theta));
	}

	results[n] = H_n;
    }

    return results;
}

void print_Complex(Complex* z)
{
    printf("(%Lf + %Lfi)", z->r, z->i);
}

// Placeholder function with similar function to 'pass' statement in Python
void pass() { }

// Calculates IFT of N provided samples.
// TODO: Understand why it doesn't work
Complex* IFT(Complex* samples, size_t N, size_t skip_n)
{
    Complex* result = (Complex*)malloc(N * sizeof(Complex));

    Complex H_n, h_k;
    long double theta;
    size_t n, k;

    for (k = 0; k < N; k++) {
	h_k.r = 0.;
	h_k.i = 0.;

	for (n = 0; n < N; n++) {
	    if (n == skip_n) {
		pass();
	    } else {
		// Calculate values
		theta = 2. * pi * n * k / N;
		H_n = samples[n];

		// (a + bi)(c + di) = (ac - bd) + (ad + bc)i
		// h_k(t_k).exp(-2.pi.n.k/N):
		h_k.r += (H_n.r * cosl(theta) - H_n.i * sinl(theta));
		h_k.i += (H_n.r * sinl(theta) + H_n.i * cosl(theta));
	    }
	}

	h_k.r = (long double)h_k.r / N;
	h_k.i = (long double)h_k.i / N;
	result[k] = h_k;

	//print_Complex(&result[k]);
	//putchar('\n');
    }

    return result;
}

// Implementation for 3b
Complex** q_3b(long double* times, size_t N)
{
    // Allocating memory for arrays to store results in
    Complex* h1_vals = (Complex*)malloc(N * sizeof(Complex));
    Complex* h2_vals = (Complex*)malloc(N * sizeof(Complex));
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));

    results[0] = h1_vals;
    results[1] = h2_vals;

    // Open file and check for success.
    FILE* fp1 = fopen("3_b_1.txt", "w");
    FILE* fp2 = fopen("3_b_2.txt", "w");
    if (!fp1 || !fp2) {
	fprintf(stderr, "Failed to open file in <q_3b>");
	exit(-1);
    }

    // Declare variables for increment, complex results, and times.
    size_t i;
    Complex y1, y2;
    long double time;

    // Print labels to file.
    fprintf(fp1, "time,h1.r,h1.i");
    fprintf(fp2, "time,h2.r,h2.i");

    // For every sample
    for (i = 0; i < N; ++i) {
	time = times[i];

	// Calculate vals of h1 & h2 at time t
	y1 = h_1(time);
	y2 = h_2(time);

	// Write to respective arrays in case we need them later...
	*(results[0] + i) = y1;
	*(results[1] + i) = y2;

	// Write y1 & y2 to file in terms of both their real and imaginary components.
	// To be plotted using Python or Gnuplot.
	fprintf(fp1, "\n%Lf, %Lf, %Lf", time, y1.r, y1.i);
	fprintf(fp2, "\n%Lf, %Lf, %Lf", time, y2.r, y2.i);
    }

    // Close files
    fclose(fp1);
    fclose(fp2);

    // Free memory
    free(times);

    return results;
}

// Questions 3.d & 3.e:
// 	- 3.d: Perform DFT on h1(t) & h2(t) --> H1(w) & H2(w)
// 	       typeof(Hn(w)) = Complex*
// 	       Calulated using DFT_functions
// 	- 3.e: Print results to screen
// Also returns array of {H1(w), H2(w)} for use in part 3.f
Complex** q_3de(long double* times, size_t N)
{
    // Allocating memory for results, H1(w), & H2(w)
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
    Complex* H1 = (Complex*)malloc(N * sizeof(Complex));
    Complex* H2 = (Complex*)malloc(N * sizeof(Complex));
}

int main()
{
    // Set number of samples & generate time values
    size_t N = 100;
    long double* times = linspace_ld(0, 2 * pi, 100);

    // Complete Q3.b, then pass the values --> array of array of Complexes
    Complex** h1_and_h2 = q_3b(times, N);

    // Compute Discrete Fourier Transforms of samples generated in Q3.b
    // DFT_samples performs operation on existing samples,
    // whereas DFT_functions generates examples based on a function passed to it.
    Complex* H1 = DFT_samples(h1_and_h2[0], times, N);
    Complex* H2 = DFT_samples(h1_and_h2[1], times, N);

    // Allocate space for & compute Inverse Fourier Transform (IFT)
    // using values from H1 & H2
    Complex* i_h1 = IFT(H1, N, 1);
    Complex* i_h2 = IFT(H2, N, 0);

    // Complex values to keep next loop readable;
    Complex H1_current;
    Complex H2_current;

    int i;
    // Printing values of DFT
    for (i = 0; i < N; i++) {
	H1_current = H1[i];
	H2_current = H2[i];

	//printf("t = %Lf\tH1_%d = %Lf + %Lfi\tH2_%d = %Lf + %Lfi\n", times[i], i, H1_current.r, H1_current.i, i, H2_current.r, H2_current.i);
    }

    FILE* fp1 = fopen("inv_1.txt", "w");
    FILE* fp2 = fopen("inv_2.txt", "w");
    if (!fp1 || !fp2) {
	fprintf(stderr, "failed opening inverse.txt");
	exit(-1);
    }

    fprintf(fp1, "time,h1.r,h1.i\n");
    fprintf(fp2, "time,h2.r,h2.i\n");
    for (i = 0; i < N; ++i) {
	fprintf(fp1, "%Lf, %Lf, %Lf\n", times[i], i_h1[i].r, i_h1[i].i);
	fprintf(fp2, "%Lf, %Lf, %Lf\n", times[i], i_h2[i].r, i_h2[i].i);
    }

    fclose(fp1);
    fclose(fp2);

    free(i_h1);
    free(i_h2);
    free(times);
    free(h1_and_h2);
    free(H1);
    free(H2);

    return 0;
}
