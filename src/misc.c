#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include <string.h>
#include <math.h>

void print_result(double *x, double *r, double *p, double *Ap, double residual, int iter, int N) {
    printf("Solution vector x*:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", x[i]);
    }
    printf("\n");

    printf("Residual vector r:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", r[i]);
    }
    printf("\n");

    printf("Search direction vector p:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", p[i]);
    }
    printf("\n");

    printf("Vector Ap:\n");
    for (int i = 0; i < N; i++) {
        printf("%f ", Ap[i]);
    }
    printf("\n");

    printf("Final residual: %f\n", residual);
    printf("Total iterations: %d\n", iter);
}

CGParams init_cg(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file\n");
        exit(EXIT_FAILURE);
    }

    CGParams params;

    if (fscanf(file, "N = %d\n", &params.N) != 1 ||
        fscanf(file, "tol = %lf\n", &params.tol) != 1 ||
        fscanf(file, "max_iter = %d\n", &params.max_iter) != 1) {
        fprintf(stderr, "Failed to read parameters\n");
        exit(EXIT_FAILURE);
    }

    params.A = (double*) malloc(params.N * params.N * sizeof(double));
    if (params.A == NULL) {
        fprintf(stderr, "Memory allocation failed for A\n");
        exit(EXIT_FAILURE);
    }

    params.b = (double*) malloc(params.N * sizeof(double));
    if (params.b == NULL) {
        fprintf(stderr, "Memory allocation failed for b\n");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "A =");
    for (int i = 0; i < params.N * params.N; i++) {
        if (fscanf(file, "%lf", &params.A[i]) != 1) {
            fprintf(stderr, "Failed to read A at index %d\n", i);
            exit(EXIT_FAILURE);
        }
    }

    fscanf(file, " b =");
    for (int i = 0; i < params.N; i++) {
        if (fscanf(file, "%lf", &params.b[i]) != 1) {
            fprintf(stderr, "Failed to read b at index %d\n", i);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
    return params;
}

void generateRandomSPDMatrix(double *A, int N) {
    int i, j;
    double sum;
    for (i = 0; i < N; i++) {
        sum = 0.0;
        for (j = 0; j < N; j++) {
            if ((rand() % 10) < 3 && i != j) {
                A[i * N + j] = A[j * N + i] = (double)(rand() % 20 + 1) - 10;
                sum += fabs(A[i * N + j]);
            } else {
                A[i * N + j] = 0.0;
            }
        }
        A[i * N + i] = sum + 1.0;
    }
}

CGParams random_init_cg(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file\n");
        exit(EXIT_FAILURE);
    }

    CGParams params;

    if (fscanf(file, "N = %d\n", &params.N) != 1 ||
        fscanf(file, "tol = %lf\n", &params.tol) != 1 ||
        fscanf(file, "max_iter = %d\n", &params.max_iter) != 1) {
        fprintf(stderr, "Failed to read parameters\n");
        exit(EXIT_FAILURE);
    }

    fclose(file);

    params.A = (double*) malloc(params.N * params.N * sizeof(double));
    if (params.A == NULL) {
        fprintf(stderr, "Memory allocation failed for A\n");
        exit(EXIT_FAILURE);
    }

    params.b = (double*) malloc(params.N * sizeof(double));
    if (params.b == NULL) {
        fprintf(stderr, "Memory allocation failed for b\n");
        exit(EXIT_FAILURE);
    }

    generateRandomSPDMatrix(params.A, params.N);
    for (int i = 0; i < params.N; i++) {
        params.b[i] = (double)(rand() % 100);
    }

    return params;
}

void print_cg_params(const CGParams params) {
    printf("Conjugate Gradient Parameters:\n");
    printf("Dimension N: %d\n", params.N);
    printf("Tolerance (tol): %g\n", params.tol);
    printf("Maximum Iterations (max_iter): %d\n", params.max_iter);

    // printf("Matrix A:\n");
    // for (int i = 0; i < params.N; i++) {
    //     for (int j = 0; j < params.N; j++) {
    //         printf("%9.4f ", params.A[i * params.N + j]);
    //     }
    //     printf("\n");
    // }

    // printf("Vector b:\n");
    // for (int i = 0; i < params.N; i++) {
    //     printf("%9.4f ", params.b[i]);
    // }
    printf("\n");
}