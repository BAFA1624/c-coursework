#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

#define SWAP(x, y)          \
    do {                    \
	typeof(x) SWAP = x; \
	x = y;              \
	y = SWAP;           \
    } while (0)

typedef struct Complex {
    double r;
    double i;
} Complex;

Complex* iniComplex(double _r, double _i)
{
    Complex* z = (Complex*)malloc(sizeof(Complex));
    z->r = _r;
    z->i = _i;
    return z;
}

double cModPtr(const Complex* z)
{
    return sqrt(z->r * z->r + z->i * z->i);
}

int compareComplex(const void* a, const void* b)
{
    const double mod_fa = cModPtr((const Complex*)a);
    const double mod_fb = cModPtr((const Complex*)b);
    if (isnan(mod_fa) || isnan(mod_fb)) {
	exit(-420);
    } else if (mod_fa > mod_fb) {
	return 1;
    } else if (mod_fa < mod_fb) {
	return -1;
    }
    return 0;
}

void cSwap(Complex* a, Complex* b)
{
    Complex tmp = *a;
    *a = *b;
    *b = tmp;
}

long partition(Complex* arr, long left_idx, long right_idx, int (*compare)(const void*, const void*))
{
    Complex x = arr[right_idx];
    long i = (left_idx - 1);
    for (long j = left_idx; j <= right_idx - 1; ++j) {
	int comparison = compare((void*)&arr[j], (void*)&x);
	// if arr[j] <= x (by complex magnitude)
	if (comparison == -1 || comparison == 0) {
	    i++;
	    cSwap(&arr[i], &arr[j]);
	}
    }
    cSwap(&arr[i + 1], &arr[right_idx]);
    return (i + 1);
}

void cQuicksort(Complex* arr, long left_idx, long right_idx, int (*compare)(const void*, const void*))
{
    // Create auxiliary stack
    long stack[right_idx - left_idx + 1];

    // initialize top of stack;
    long top = -1;

    // Push initial values to stack
    stack[++top] = left_idx;
    stack[++top] = right_idx;

    // Keep popping from stack while it isn't empty
    while (top >= 0) {
	// Pop right_idx & left_idx
	right_idx = stack[top--];
	left_idx = stack[top--];

	// Set pivot element at correct position in sorted array
	long pivot = partition(arr, left_idx, right_idx, compare);

	// If elements on left side of pivot, push left side to stack
	if (pivot - 1 > left_idx) {
	    stack[++top] = left_idx;
	    stack[++top] = pivot - 1;
	}

	// If elements on right side of pivot, push right side to stack
	if (pivot + 1 < right_idx) {
	    stack[++top] = pivot + 1;
	    stack[++top] = right_idx;
	}
    }
}

int main()
{
}
