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
long double c_mod(Complex z)
{
    return sqrtl(z.r * z.r + z.i * z.i);
}
long double c_mod_ptr(Complex* z)
{
    return c_mod(*z);
}

// Calculate argument of complex number
long double c_arg(Complex z)
{
    return atan2l(z.i, z.r);
}
long double c_arg_ptr(Complex* z)
{
    return atan2l(z->i, z->r);
}

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

// TODO: Implement IFT once confirmed output of DFT
Complex* IFT(Complex* samples, int N)
{
    Complex* results = (Complex*)malloc(N * sizeof(Complex));

    Complex H_n, h_k;
    long double theta;
    int n, k;

    for (k = 0; k < N; k++) {
	h_k.r = 0.;
	h_k.i = 0.;

	for (n = 0; n < N; n++) {

	    // Calculate values
	    theta = 2. * pi * n * k / N;
	    H_n = samples[k];

	    // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
	    // h_k(t_k).exp(-2.pi.n.k/N):
	    h_k.r += (H_n.r * cosl(theta) - H_n.i * sinl(theta));
	    h_k.i += (H_n.r * sinl(theta) + H_n.i * cosl(theta));
	}

	h_k.r = (long double)h_k.r / N;
	h_k.i = (long double)h_k.i / N;
	printf("h_k.r = %Lf, h_k.i = %Lf\n", h_k.r, h_k.i);
	results[n] = h_k;
    }

    return results;
}

// Implementation for 3b
Complex** q_3b(int N)
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
    int i;
    Complex y1, y2;
    long double* times = linspace_ld(0., 2 * pi, N);
    long double time;

    // Print labels to file.
    fprintf(fp1, "time,h1.r,h1.i");
    fprintf(fp2, "time,h2.r,h2.i");

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

int main()
{
    int N = 100;

    Complex** h1_h2 = q_3b(N);

    long double* times = linspace_ld(0, 2 * pi, N);
    Complex* H1 = DFT_functions(h_1, times, N);
    Complex* H2 = DFT_functions(h_2, times, N);

    Complex* i_h1 = IFT(H1, 100);
    Complex* i_h2 = IFT(H2, 100);

    Complex H1_current;
    Complex H2_current;

    int i;
    for (i = 0; i < N; i++) {
	H1_current = H1[i];
	H2_current = H2[i];

	printf("t = %Lf\tH1_%d = %Lf + %Lfi\tH2_%d = %Lf + %Lfi\n", times[i], i, H1_current.r, H1_current.i, i, H2_current.r, H2_current.i);
    }

    FILE* fp = fopen("inverse.txt", "w");
    if (!fp) {
	fprintf(stderr, "failed opening inverse.txt");
	exit(-1);
    }

    fprintf(fp, "time,h1.r,h1.i\n");
    for (i = 0; i < N; ++i) {
	fprintf(fp, "%Lf, %Lf, %Lf\n", times[i], i_h1[i].r, i_h1[i].i);
    }

    fclose(fp);

    free(i_h1);
    free(i_h2);
    free(times);
    free(h1_h2);
    free(H1);
    free(H2);

    return 0;
}
