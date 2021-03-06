#include <stdio.h>
#include <math.h>
#include <stdlib.h>

long double pi = 3.14159265358979323846;
long double pi_2 = 1.57079632679489661923;
long double pi_4 = 0.78539816339744830962;
long double e = 2.7182818284590452354;


// Implementing some simple complex number functionality

typedef struct Complex 
{    
    long double r;      // Real part;
    long double i;      // Imaginary part;
} Complex;
// Calculate square of modulus for complex numbers.
long double c_mod(Complex z) 
{
    return sqrtl(z.r * z.r + z.i * z.i);
}
long double c_mod_ptr(Complex *z)
{
    return c_mod(*z);
}

// Version of c_arg using atan2l

long double c_arg(Complex z)
{
    return atan2l(z.i, z.r);
}
long double c_arg_ptr(Complex *z)
{
    return atan2l(z->i, z->r);
}

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

void q_3b(int N)
{
    // Open file and check for success.
    FILE *fp1 = fopen("3_b_1.txt", "w");
    FILE *fp2 = fopen("3_b_2.txt", "w");
    if (!fp1 || !fp2)
    {
        fprintf(stderr, "Failed to open file in <q_3b>");
        exit(-1);
    }

    // Decl + initialize vars for time increment, time interval,
    // number of vars N, Complex numbers y1 & y2, increment i
    long double t_ini = 0.;
    int i;
    long double dt = (long double)2 * pi / N;
    Complex y1, y2;

    // Print labels to file.
    fprintf(fp1, "time,h1.r,h1.i");
    fprintf(fp2, "time,h2.r,h2.i");

    for (i = 0; i < N; ++i)
    {
	// Calculate vals of h1 & h2 at time t
        y1 = h_1(t_ini);
        y2 = h_2(t_ini);

	// Write y1 & y2 to file in terms of both their real and imaginary components.
	// To be plotted using Python or Gnuplot.
        fprintf(fp1, "\n%Lf, %Lf, %Lf", t_ini, y1.r, y1.i);
	fprintf(fp2, "\n%Lf, %Lf, %Lf", t_ini, y2.r, y2.i);  
        
	// Increment time
	t_ini = t_ini + dt;
    }

    // Close file
    fclose(fp1);
    fclose(fp2);
}

int main()
{
    q_3b(100);
    return 0;
}
