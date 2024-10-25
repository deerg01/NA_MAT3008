#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "nr.h"  
/*
cc -I ./NRs/ansi/other -o run ./Homework#4/Code/test.c ./NRs/ansi/recipes/ran1.c ./NRs/ansi/recipes/gasdev.c
*/
void generate_uniform(int sample, const char *filename, double a, double b, long *idum) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int j = 0; j < sample; ++j) {
        fprintf(file, "%lf\n", ran1(idum) * (b - a) + a); 
    }

    fclose(file);
}
void generate_gaussian(int sample, const char *filename, double m, double s, long *idum) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int j = 0; j < sample; ++j) {
        fprintf(file, "%lf\n", gasdev(idum) * s + m); 
    }
    
    fclose(file);
}
int main(void) {
    long idum = time(0);
    
    int n[] = {100, 1000, 10000, 100000};
    int n_samples = sizeof(n) / sizeof(n[0]);
    
    double a = -3, b = 4, m = 0.5, s = 1.5;

    for (int i = 0; i < n_samples; ++i) {
        int sample = n[i];

        char uf[15], gf[15];
        sprintf(uf, "u_%d.txt", sample);
        sprintf(gf, "g_%d.txt", sample);

        generate_uniform(sample, uf, a, b, &idum);
        generate_gaussian(sample, gf, m, s, &idum);
    }

    return 0;
}
