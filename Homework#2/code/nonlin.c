#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define XACC 1.0e-6
/*
cc -o run2 ./Homework#2/code/nonlin.c ./NRs/ansi/recipes/rtbis.c ./NRs/ansi/recipes/rtflsp.c ./NRs/ansi/recipes/rtsec.c ./NRs/ansi/recipes/rtnewt.c ./NRs/ansi/recipes/rtsafe.c ./NRs/ansi/recipes/zbrak.c ./NRs/ansi/other/nrutil.c
*/

//copying from NR source files
float rtbis(float (*func)(float), float x1, float x2, float xacc);
float rtflsp(float (*func)(float), float x1, float x2, float xacc);
float rtsec(float (*func)(float), float x1, float x2, float xacc);
float rtnewt(void (*funcd)(float, float *, float *), float x1, float x2, float xacc);
float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2, float xacc);
void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[], float xb2[], int *nb);

//written by me
float sgn(float value) { //sign function
    return (value > 0) - (value < 0);  
}
float muller(float (*func)(float), float x1, float x2, float xacc) {
    float p1 = x1, p2 = x2, p3;
    float p0 = (p1 + p2) / 2; //median
    float f0, f1, f2, h1, h2, a, b, c, rad, den;
    int iter, max_iter = 100; 
    
    for (iter = 0; iter < max_iter; iter++) {
        f0 = func(p0);  
        f1 = func(p1);  
        f2 = func(p2);  
        //get a,b,c
        c = f2;
        b = ((p0 - p2) * (p0 - p2) * (f1 - f2) - (p1 - p2) * (p1 - p2) * (f0 - f2)) /
            ((p0 - p2) * (p1 - p2) * (p0 - p1));
        a = ((p1 - p2) * (f0 - f2) - (p0 - p2) * (f1 - f2)) /
            ((p0 - p2) * (p1 - p2) * (p0 - p1));
        //Discriminant, root formulas
        rad = sqrt(b * b - 4 * a * c);
        den = b + sgn(b) * rad;
        p3 = p2 - (2 * c) / den;
        //end of iteration
        if (fabs(p3 - p2) < xacc) {
            return p3;
        }
        p0 = p1;
        p1 = p2;
        p2 = p3;
    }

    //error and escape
    printf("Muller method failed to converge\n");
    return p3;
}
float nonlin(float x){
    return sin(x) - 0.5 * x + cos(5 * x);
}
void nonlin_(float x, float *f, float *df){ //derivative
    *f = sin(x) - 0.5 * x + cos(5 * x);
    *df = cos(x) - 0.5 - 5 * sin(5 * x);
}

int main(int argc, char const *argv[]){
    float xb1[101], xb2[101]; 
    int i, nb = 20; 
    float x1 = 1.0, x2 = 10.0; //interval [1.0, 10.0]

    zbrak(nonlin, x1, x2, 100, xb1, xb2, &nb);

    for (i = 1; i <= nb; i++) { 
        x1 = xb1[i];
        x2 = xb2[i];
        printf("*** Range [%f, %f] ***\n\n", x1, x2); 
        printf("Bisection method: %f\n\n", rtbis(nonlin, x1, x2, XACC));
        printf("Linear Interpolation method: %f\n\n", rtflsp(nonlin, x1, x2, XACC));
        printf("Secant method: %f\n\n", rtsec(nonlin, x1, x2, XACC));
        printf("Newton-Raphson method: %f\n\n", rtnewt(nonlin_, x1, x2, XACC)); 
        printf("Newton with bracketing: %f\n\n", rtsafe(nonlin_, x1, x2, XACC)); 
        printf("Muller Method: %f\n\n", muller(nonlin, x1, x2, XACC));
    }
    printf("********************************\n");
}