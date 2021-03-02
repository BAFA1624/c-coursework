#include <stdio.h>
#include <math.h>

#define M_PI		3.14159265358979323846
#define M_PI_2		1.57079632679489661923
#define M_PI_4		0.78539816339744830962


// Implementing some simple complex number functionality

typedef struct Complex 
{    
    long double r;      // Real part;
    long double i;      // Imaginary part;
} Complex;
// Calculate square of modulus for complex numbers.
long double c_mod_sq(Complex z) 
{
    return z.r * z.r + z.i * z.i;
}
long double c_mod_sq_ptr(Complex *z)
{
    return c_mod_sq(*z);
}
/*
// Version of c_arg using atan2l

long double c_arg(Complex z)
{
    return atan2l(z.i, z.r);
}
long double c_arg_ptr(Compmlex *z)
{
    return atan2l(z->i, z->r);
}

*/
// atanl is long double arc tan.
// Calculate argument for complex numbers.
long double c_arg(Complex z)
{    
    switch (z.r > 0)
    {
        // Case 1 [true]: real part > 0
        case 1: 
            return atanl(z.i / z.r);
        // Case 0 [false]: real part <= 0
        case 0:
            if (z.r != 0)
            {
                // Use ternary operator to return correct expression
                return (z.i < 0) ? atanl(z.i / z.r) - M_PI : atanl(z.i / z.r) + M_PI;
            } 
            else
            {
                if (z.i == 0)
                {
                    fprintf(stderr, "Undefined value in <c_arg>, z.r = %ld, z.i = %ld.", z.r, z.i);
                    exit(-1);
                }
                else
                {
                    return (z.i > 0) ? M_PI_2 : -M_PI_2;
                }
            }
        default:
            printf("\n[Defaulted in c_arg]\n");
    }
}
long double c_arg_ptr(Complex *z)
{
    return c_arg(*z);
}


Complex h_1()
