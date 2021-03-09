#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

long double pi = 3.14159265358979323846;
long double pi_2 = 1.57079632679489661923;
long double pi_4 = 0.78539816339744830962;
long double e = 2.7182818284590452354;

// Declarations for functions lacking them in header files

#define errno_t int
#define rsize_t size_t
errno_t memcpy_s(void* restrict dest, rsize_t destsz, const void* restrict src, rsize_t count);

// -----------------------------------------------------
// Implementing some simple complex number functionality

typedef struct Complex {
    long double r; // Real part;
    long double i; // Imaginary part;
} Complex;

// Calculate square of modulus for complex numbers.
long double c_mod(const Complex z) { return sqrtl(z.r * z.r + z.i * z.i); }
long double c_mod_ptr(const Complex* z) { return c_mod(*z); }

// Calculate argument of complex number
long double c_arg(const Complex z) { return atan2l(z.i, z.r); }
long double c_arg_ptr(const Complex* z) { return atan2l(z->i, z->r); }

void print_Complex(const Complex* z)
{
    printf("(%Lf + %Lfi)", z->r, z->i);
}

// End of complex number functions
// ------------------------------------------------------
// Simple struct for returning data in 3j+

typedef struct Measurement {
    int* n_arr;
    long double* times;
    Complex* z_arr;
} Measurement;

// End of Measurement struct/functions
// -------------------------------------------------------
// General functionality:
// Functions which facilitate completing the questions

// Generic write error message, then exit execution
void error_exit(const char* err_msg)
{
    fprintf(stderr, err_msg);
    exit(-1);
}

// Implementation of h1 & h2
Complex h_1(const long double time)
{
    // Calculate h1 using Euler's formula.
    Complex result = {
	.r = cosl(time) + cosl(5 * time),
	.i = sinl(time) + sinl(5 * time)
    };
    return result;
}
Complex h_2(const long double time)
{
    // calculate exponent, theta, of h2
    long double theta = (time - pi) * (time - pi) / 2;
    // Initialize Complex number with result.
    Complex result = {
	.r = powl(e, theta),
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
	error_exit("\nrealloc failed in arr_pop_ld\n");
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
	error_exit("realloc failed in arr_pop_Complex\n");
    }
    return arr;
}

// Produces range of N long doubles from start -> end exclusive. Simple version of numpy.linspace
long double* linspace_ld(const long double start, const long double end, const size_t N)
{
    // Allocate memory for result, calculate increment size for N values in given range
    long double* arr = (long double*)malloc(N * sizeof(long double));
    if (!arr) {
	error_exit("\nmalloc failed in linspace_ld\n");
    }

    // Calculate increment for N evenly spaced values from start -> end
    // 'val' keeps track of current value
    long double increment = (double)(end - start) / N;
    long double val = start;

    // Add values to arr, incrementing val every time
    size_t i;
    for (i = 0; i < N; i++) {
	arr[i] = val;
	val += increment;
    }

    return arr;
}

// Produces range of N complex values for given function func for provided times
Complex* linspace_Complex(Complex (*func)(long double), const long double* times, const size_t N)
{
    // Allocate memory for result, error if fails
    Complex* arr = (Complex*)malloc(N * sizeof(Complex));
    if (!arr) {
	error_exit("\nmalloc failed in linspace_Complex\n");
    }

    // Generate a value in arr for every provided value in times
    size_t i;
    for (i = 0; i < N; ++i) {
	arr[i] = func(times[i]);
    }
    return arr;
}

// Discrete Fourier Transform taking array of samples
Complex* DFT_samples(const Complex* samples, const size_t N)
{
    // Alloc memory for resulting array & declare vars
    Complex* arr = (Complex*)malloc(N * sizeof(Complex));
    if (!arr) {
	error_exit("\nmalloc failed in DFT_samples\n");
    }

    // H_n & h_k keep track of which values are currently in use
    // theta stores value of exponent h_k(t_k).exp(-2.pi.n.k/N) for use in Euler's formula
    Complex H_n, h_k;
    long double theta;
    size_t n, k;

    // For every values H_n, sum contributions of all h_k then add to resulting array
    for (n = 0; n < N; n++) {
	// initialize H_n
	H_n.r = 0.;
	H_n.i = 0.;

	for (k = 0; k < N; k++) {
	    // Calculate exponent of h_k(t_k).exp(-2.pi.n.k/N)
	    theta = 2. * pi * n * k / N;
	    h_k = samples[k];

	    // (a + bi)(c + di) = (ac - bd) + (ad + bc)i
	    // h_k(t_k).exp(-2.pi.n.k/N):
	    H_n.r += (h_k.r * cosl(theta) + h_k.i * sinl(theta));
	    H_n.i += (h_k.i * cosl(theta) - h_k.r * sinl(theta));
	}

	arr[n] = H_n;
    }

    return arr;
}

// Discrete Fourier Transform taking sample function
Complex* DFT_function(Complex (*func)(long double), const long double* times, const int N)
{
    // Takes:
    // 	- func: Ptr to function to perform DFT on.
    // 	- times: Long double array of sample times.
    // 	- N: Size of times, incidentally also number of samples.
    // Does:
    //  - Generates values for func using time values provided in times.
    //  - Passes generated values --> DFT_samples.
    //  - Returns result

    Complex* samples = linspace_Complex(func, times, N);
    Complex* arr = DFT_samples(samples, N);

    return arr;
}

// Placeholder function with similar function to 'pass' statement in Python
// Does nothing by design, e.g. for use in if/else or switch statements
void pass() { }

// Calculates IFT of N provided samples.
Complex* IFT(Complex* samples, size_t N, size_t skip_n)
{
    // Alloc memory for resulting array & declare vars
    Complex* arr = (Complex*)malloc(N * sizeof(Complex));
    if (!arr) {
	error_exit("\nmalloc failed in IFT\n");
    }

    // H_n & h_k keep track of which values are currently in use
    // theta stores value of exponent h_k(t_k).exp(-2.pi.n.k/N) for use in Euler's formula
    Complex H_n, h_k;
    long double theta;
    size_t n, k;

    // For every values H_n, sum contributions of all h_k then add to resulting array
    for (k = 0; k < N; k++) {
	// initialize h_k
	h_k.r = 0.;
	h_k.i = 0.;

	for (n = 0; n < N; n++) {
	    if (n == skip_n) {
		pass();
	    } else {
		// Calculate exponent of H_n(W_n).exp(-2.pi.n.k/N)
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

	arr[k] = h_k;
    }

    return arr;
}

// Compare function for sorting values in 3.k using qsort from stdlib.h
// Takes const (void*) cast ptrs to two Complex numbers
// qsort(void* base, size_t nitems, size_t size, int (*compar)(const vooid*, const void*))
// int return tells qsort which value is >, ==, or < the other
// result > 0 --> a > b
// result < 0 --> b > a
// result = 0 --> a = b
int compare_Complex(const void* a, const void* b)
{
    // Calculate amplitude of Complex vars
    const long double mod_fa = c_mod_ptr((const Complex*)a);
    const long double mod_fb = c_mod_ptr((const Complex*)b);

    // Comparison of two Complex numbers
    // If nan values detected --> error out; else --> return appropriate value for comparison
    if (isnan(mod_fa) || isnan(mod_fb)) {
	    error_exit("\nnan values in compare_Complex.\n");
    } else if (mod_fa > mod_fb) {
	    return 1;
    } else if (mod_fa < mod_fb) {
	    return -1;
    } else {
	    return 0;
    }

    return 0;
}

// End of general functions
// -----------------------------------------------------
// Question answers

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

	// Write to respective arrays to they can be used in the next questions
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
    size_t i;
    // Allocating memory for results; will store {H1, H2}
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
    if (!results) {
	error_exit("\nmalloc failed allocating results in q_3de\n");
    }

    // Calculating H1 & H2 using DFT_function
    Complex* H1 = DFT_function(h_1, times, N);
    Complex* H2 = DFT_function(h_2, times, N);

    // Add H1 & H2 --> results
    results[0] = H1;
    results[1] = H2;

    // Printing H1 & H2 using custom print_Complex function
    printf("\nPrinting H1:\n");
    for (i = 0; i < N; ++i) {
	printf("t = %Lf\tH1 = ", times[i]);
	print_Complex(&H1[i]);
	putchar('\n');
    }
    printf("\nPrinting H2:\n");
    for (i = 0; i < N; ++i) {
	printf("t = %LF\tH2 = ", times[i]);
	print_Complex(&H2[i]);
	putchar('\n');
    }

    return results;
}

// Applies Inverse Fourier Transform (IFT) to provided sets of H_x in samples.
// In this case samples will provided by result of q_3de.
Complex** q_3f(Complex** samples, size_t N)
{
    // Unpack H1 & H2 from samples.
    Complex* H1 = samples[0];
    Complex* H2 = samples[1];

    // Calculate Inverse Fourier Transform of H1 & H2 respectively.
    Complex* h1_prime = IFT(H1, N, 1);
    Complex* h2_prime = IFT(H2, N, 0);

    // Allocate memory for results, check for fail.
    Complex** results = (Complex**)malloc(2 * sizeof(Complex*));
    if (!results) {
	error_exit("\nFailed malloc for results in q_3f.\n");
    }

    // Pack h1_prime & h2_prime into results so it can be returned for q_3f
    results[0] = h1_prime;
    results[1] = h2_prime;

    return results;
}

void q_3g(Complex** data, long double* times, size_t N)
{
    // Open files to write data to. Check for file open error.
    FILE* fp1 = fopen("inv_1.txt", "w");
    FILE* fp2 = fopen("inv_2.txt", "w");
    if (!fp1 || !fp2) {
	error_exit("\nFailed opening file in q_3g.\n");
    }

    // Unpack data
    Complex* h1_prime = data[0];
    Complex* h2_prime = data[1];

    // Print file headers to inv_1.txt & inv_2.txt
    fprintf(fp1, "time,real,imag\n");
    fprintf(fp2, "time,real,imag\n");

    // Write data to file
    size_t i;
    for (i = 0; i < N; ++i) {
	fprintf(fp1, "%Lf,%Lf,%Lf\n", times[i], h1_prime[i].r, h1_prime[i].i);
	fprintf(fp2, "%Lf,%Lf,%Lf\n", times[i], h2_prime[i].r, h2_prime[i].i);
    }

    // Closing files
    fclose(fp1);
    fclose(fp2);
}

// Read data in from data.txt, return as array of Complexes
Measurement q_3h(const char* filename, size_t N)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) {
	    error_exit("\nFailed to open data.txt in q_3h.\n");
    }

    // Allocating vars to keep track of array size and arrays for data
    size_t sz = 0;
    long double* times = (long double*)malloc(N * sizeof(long double));
    if (!times) {
	    error_exit("\nmalloc failed for times in q_3h.\n");
    }

    int* n_arr = (int*)malloc(N * sizeof(int));
    if (!times) {
	    error_exit("\nmalloc failed for n_arr in q_3h.\n");
    }

    Complex* z_arr = (Complex*)malloc(N * sizeof(Complex));
    if (!z_arr) {
	    error_exit("\nmalloc failed for z_arr in q_3h.\n");
    }

    long double real = 0.;
    long double imag = 0.;

    // Use fscanf to read from file
    while (fscanf(fp, "%d, %Lf, %Lf, %Lf", &n_arr[sz], &times[sz], &real, &imag) != EOF && sz < N) {
	    z_arr[sz].r = real;
	    z_arr[sz].i = imag;
	    sz++;
    }

    fclose(fp);

    // Create result & pack it with data from file read
    Measurement result;
    result.n_arr = n_arr;
    result.times = times;
    result.z_arr = z_arr;

    return result;
}

// Apply DFT to h3
Complex* q_3j(const Complex* data, const size_t N)
{
    // Applying DFT to data with DFT_samples
    Complex* result = DFT_samples(data, N);
    return result;
}



// Applies IFT to 4 terms of H3 with largest amplitude
/*Complex* q_3k(const Complex* samples, const size_t N)
{
    // Allocating memory for indexes of terms with the four largest amplitudes
    size_t* idx = (size_t*)malloc(4 * sizeof(size_t));
    size_t i, max_idx = 0;

    // Copying samples to array for sorting
    Complex* sorted_samples = (Complex*)malloc(N * sizeof(Complex));
    errno_t status = memcpy_s((void*)sorted_samples, N * sizeof(Complex), samples, N);
    if (!status) {
	error_exit("\nmemcpy_s in q_3k returned failing status.\n");
    }

    // Sorting samples
    qsort((void*)sorted_samples, N, sizeof(Complex), compare_Complex);

    //
    
}*/

// End of question answers
// -----------------------------------------------------
// main()

int main()
{
    // Set number of samples & generate time values
    size_t i, N = 100;
    long double* times = linspace_ld(0, 2 * pi, N);

    // Complete Q3.b, then pass the values --> array of array of Complexes
    Complex** h1_and_h2 = q_3b(times, N);

    // Complete Q3.d & Q3.e, then pass the values --> array of array of Complexes
    Complex** H1_and_H2 = q_3de(times, N);

    Complex** prime_h1_and_h2 = q_3f(H1_and_H2, N);

    q_3g(prime_h1_and_h2, times, N);

    // Data provided uses N=200
    N = 200;

    Measurement h3 = q_3h("data.txt", N);

    Complex* H3 = q_3j(h3.z_arr, N);

    for (i = 0; i < N; ++i) {
        printf("\nH3_%ld = ", i);
        print_Complex(&H3[i]);
        putchar('\n');
    }

    /*
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
    */

    return 0;
}
