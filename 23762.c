#include <stdio.h>
#include <math.h>
#include <stdlib.h>

long double pi = 3.1415926;
long double pi_2 = 1.5707963;
long double pi_4 = 0.7853981;
long double e = 2.7182818;


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
                return (z.i < 0) ? atanl(z.i / z.r) - pi : atanl(z.i / z.r) + pi;
            } 
            else
            {
                if (z.i == 0)
                {
                    fprintf(stderr, "Undefined value in <c_arg>, z.r = %Lf, z.i = %Lf.", z.r, z.i);
                    exit(-1);
                }
                else
                {
                    return (z.i > 0) ? pi_2 : -pi_2;
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


Complex h_1(long double time)
{
    Complex result = {
        .r = cosl(time) + cosl(5 * time),
        .i = sinl(time) + sinl(5 * time)
    };
    return result; 
}
Complex h_2(long double time)
{
    long double x = (time - pi) * (time - pi) / 2;
    Complex result = {
        .r = powl(e, x),
        .i = 0
    };
    return result;
}

void q_3b(char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        fprintf(stderr, "Failed to open file \"%s\" in <q_3b>", filename);
        exit(-1);
    }
    long double t_ini = 0.;
    int i, N = 100;
    long double dt = (long double)2 * pi / N;
    printf("\nt_ini = %Lf, dt = %Lf\n", t_ini, dt);
    Complex y1, y2;
    for (i = 0; i < N; ++i)
    {
        y1 = h_1(t_ini);
        y2 = h_2(t_ini);
        // printf("\n|h_1(%ld)| = |%ld %ld| = %ld, |h_2(%ld)| = |%ld %ld| = %ld\n", t_ini, y1.r, y1.i, c_mod(y1), t_ini, y2.r, y2.i, c_mod(y2));
        fprintf(fp, "%Lf, %Lf, %Lf\n", t_ini, c_mod(y1), c_mod(y2));
        t_ini = t_ini + dt;
    }
    fclose(fp);
    putchar('\n');
}

int main()
{
    printf("\npi = %Lf\n", e);
    q_3b("test.txt");
    return 0;
}
