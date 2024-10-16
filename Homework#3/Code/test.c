#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nr.h"
#include "nrutil.h"

/*
cc -I ./NRs/ansi/other -o run ./Homework#3/Code/test.c ./NRs/ansi/recipes/ludcmp.c ./NRs/ansi/recipes/gaussj.c ./NRs/ansi/recipes/svdcmp.c ./NRs/ansi/other/nrutil.c ./NRs/ansi/recipes/lubksb.c ./NRs/ansi/recipes/pythag.c ./NRs/ansi/recipes/svbksb.c ./NRs/ansi/recipes/mprove.c
*/
void read(const char *fn, float ***A, float **b, int *n) {
    FILE *file = fopen(fn, "r");
    if (file == NULL) {
        perror("error: cannot open file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d", n);
    fscanf(file, "%d", n);
    

    //matrix A
    *A = (float **)malloc(((*n)+1) * sizeof(float *));
    for (int i = 1; i < *n+1; i++) {
        (*A)[i] = (float *)malloc(((*n)+1) * sizeof(float)); 
    }
    for (int i = 1; i < ((*n)+1); i++) {
        for (int j = 1; j < ((*n)+1); j++) {
            fscanf(file, "%f", &((*A)[i][j]));  
        }
    }
    //vector b
    *b = (float *)malloc(((*n)+1) * sizeof(float));
    for (int i = 1; i < *n+1; i++) {
        fscanf(file, "%f", &((*b)[i])); 
    }
    fclose(file);
}
void read_check(float **A, float *b, int n){
    printf("%d\n", n);
    printf("Matrix A:\n");
    for (int i = 1; i < n+1; i++) {
        for (int j = 1; j < n+1; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }

    printf("Vector b:\n");
    for (int i = 1; i < n+1; i++) {
        printf("%f ", b[i]);
    }
    printf("\n");
}
void gj(float **A, float *b, int n){ //Gauss-Jordan elimination
    float **A_ = (float **)malloc((n + 1) * sizeof(float *));
    for (int i = 1; i <= n; i++) {
        A_[i] = (float *)malloc((n + 1) * sizeof(float));
        memcpy(A_[i], A[i], (n + 1) * sizeof(float));
    }
    //2차원 배열의 b가 필요
    float **b_ = (float**)malloc((n+1) * sizeof(float*));
    for (int i = 0; i < n+1; i++) {
        b_[i] = (float*)malloc((n+1) * sizeof(float));
        for (int j = 0; j < n+1; j++){
            b_[i][j] = 0;
        } 
        b_[i][1] = b[i]; 
    }
    b_[0][1] = 0; 

    gaussj(A_, n, b_, 1);

    //printing output
    //for (int i = 1; i < n+1; i++){
    //    printf("%f ", b[i]);
    //} printf("\n");
    for (int i = 1; i < n+1; i++){
        printf("%f ", b_[i][1]);
    } printf("\n");

    //deallocate
    for (int i = 0; i < n+1; i++) {
        free(b_[i]);
    }
    free(b_); 
}
void lu(float **A, float *b, int n){ //LU decomposition
    int *idx = ivector(1, n);
    float d; 

    float **A_ = (float **)malloc((n + 1) * sizeof(float *));
    for (int i = 1; i <= n; i++) {
        A_[i] = (float *)malloc((n + 1) * sizeof(float));
        memcpy(A_[i], A[i], (n + 1) * sizeof(float));
    }
    float *b_ = (float *)malloc((n + 1) * sizeof(float)); 
    for (int i = 1; i <= n; i++) {
        b_[i] = b[i]; 
    }

    ludcmp(A_, n, idx, &d);
    lubksb(A_, n, idx, b_);

    for (int i = 1; i <= n; i++) {
        printf("%f ", b_[i]);
    }
    printf("\n");

    //deallocate 
    for (int i = 1; i <= n; i++) {
        free(A_[i]);
    }
    free(A_);
    free(b_);
    free_ivector(idx, 1, n); 
}
void lu_(float **A, float *b, int n){ //LU decomposition with iterative improvement
    int *idx = ivector(1, n);
    float d;

    float **A_ = (float **)malloc((n + 1) * sizeof(float *));
    for (int i = 1; i <= n; i++) {
        A_[i] = (float *)malloc((n + 1) * sizeof(float));
        memcpy(A_[i], A[i], (n + 1) * sizeof(float));
    }
    float *b_ = (float *)malloc((n + 1) * sizeof(float));
    memcpy(b_, b, (n + 1) * sizeof(float));

    ludcmp(A_, n, idx, &d);
    lubksb(A_, n, idx, b_);

    mprove(A, A_, n, idx, b, b_);
    
    for (int i = 1; i <= n; i++) {
        printf("%f ", b_[i]);
    }
    printf("\n");

    // deallocate
    for (int i = 1; i <= n+1; i++) {
        free(A_[i]);
    }
    free(A_);
    free(b_);
    free_ivector(idx, 1, n);
}
void svd(float **A, float *b, int n){ //Singular Value Decomposition
    float *w = vector(1, n);     // Singular values(sigma)
    float **v = matrix(1, n, 1, n); // V matrix
    float **A_ = (float **)malloc((n + 1) * sizeof(float *));
    for (int i = 1; i <= n; i++) {
        A_[i] = (float *)malloc((n + 1) * sizeof(float));
        memcpy(A_[i], A[i], (n + 1) * sizeof(float));
    }

    svdcmp(A_, n, n, w, v);
    
    for (int i = 1; i <= n; i++) {
        if (w[i] < 1e-6) w[i] = 0;
    }

    float *x = vector(1, n);  //solution vector

    svbksb(A_, w, v, n, n, b, x); //get

    //printing output
    for (int i = 1; i <= n; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");

    // deallocate
    free_vector(w, 1, n);
    free_vector(x, 1, n);
    free_matrix(v, 1, n, 1, n);
}

void detNInverse(float **A, int n){
    float d; //determinant
    float **Ai = matrix(1, n, 1, n); //inverse matrix
    int *idx = ivector(1, n) ;
    float **A_ = (float **)malloc((n + 1) * sizeof(float *));
    for (int i = 1; i <= n; i++) {
        A_[i] = (float *)malloc((n + 1) * sizeof(float));
        memcpy(A_[i], A[i], (n + 1) * sizeof(float));
    }

    ludcmp(A_, n, idx, &d);
    //get determinant
    for (int i = 1; i <= n; i++) {
        d *= A_[i][i];
        //printf("determinant: %f\n", d);
    }
    printf("determinant: %f\n", d);

    //get inverse matrix
    for (int i = 1; i <= n; i++) {
        float *col = vector(1, n);
        for (int j = 1; j <= n; j++) col[j] = 0.0;
        col[i] = 1.0;
        lubksb(A_, n, idx, col);
        for (int j = 1; j <= n; j++) Ai[j][i] = col[j];
        free_vector(col, 1, n);
    }
    printf("Inverse matrix:\n");
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            printf("%f ", Ai[i][j]);
        }
        printf("\n");
    }
}
int main() {
    char name[20];
    float **A, *b;
        int n = 0;
    for(int i = 1; i < 4; i++){
        sprintf(name, "lineq%d.dat", i);
        printf("***** result for %s *****\n\n", name);
        //data reading
        read(name, &A, &b, &n);
        //read_check(A, b, n);
    
        //task1
        //printf("%f %f %d\n", A[0][0], b[0], n);
        printf("Gauss-Jordan Elimination: \n");
        gj(A, b, n);
        printf("LU decomposition: \n");
        lu(A, b, n);
        printf("Singular Value Decomposition: \n");
        svd(A, b, n);
        printf("\n");

        //task2
        printf("LU decomposition with mprove(): \n");
        lu_(A, b, n);
        printf("\n");

        //task3
        detNInverse(A, n);
        printf("\n");

    }
    //memory deallocation
        for (int i = 0; i < n; i++) {
            free(A[i]);
        }
        free(A);
        free(b);

    return 0;
}