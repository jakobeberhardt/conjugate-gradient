#include <stdio.h>
#include <stdlib.h>
#include "misc.h"
#include <string.h>
#include <math.h>
#include <time.h>

#define SEED 1234
#define MAX_PRINT 4

int getRowCount(int totalRows, int size, int rank) {
	return (totalRows / size) + (totalRows % size > rank);
}

double wtime() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + ts.tv_nsec / 1E9;
}

void print_result(double *x, double *r, double *p, double *Ap, double residual, int iter, int N, int conv) {
    if(conv) {
        printf("Found sufficient solution in iteration %d with residual of %.13f\n", iter, residual);
        printf("\n");
        printf("Solution vector x*:\n");
        for (int i = 0; i < MAX_PRINT; i++) {
            printf("%f ", x[i]);
        }

    } else {
        printf("Stopped in iteration %d with residual of %.13f\n", iter, residual);
        printf("\n");
        printf("Solution vector x%d:\n", iter);
        for (int i = 0; i < MAX_PRINT; i++) {
            printf("%f ", x[i]);
        }
    }
    printf(" ...");
    printf("\n");
    printf("\n");
    printf("Residual vector r:\n");
    for (int i = 0; i < MAX_PRINT; i++) {
        printf("%f ", r[i]);
    }
    printf(" ...");
    printf("\n");
    printf("\n");

    printf("Search direction vector p:\n");
    for (int i = 0; i < MAX_PRINT; i++) {
        printf("%f ", p[i]);
    }
    printf(" ...");
    printf("\n");

    // printf("Vector Ap:\n");
    // for (int i = 0; i < MAX_PRINT; i++) {
    //     printf("%f ", Ap[i]);
    // }
    // printf("\n");
}

CGParams read_cg_config(const char* filename) {
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
    return params;
}

void load_problem(const char* filename, double *A, double *b, CGParams params) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        fprintf(stderr, "Error opening file\n");
        exit(EXIT_FAILURE);
    }   

    if (A == NULL) {
        fprintf(stderr, "Memory allocation failed for A\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    if (b == NULL) {
        fprintf(stderr, "Memory allocation failed for b\n");
        free(A);
        fclose(file); 
        exit(EXIT_FAILURE);
    }

    if (fscanf(file, "A =") != 0) {
        fprintf(stderr, "Failed to read 'A ='\n");
        free(A);
        free(b);
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < params.N * params.N; i++) {
        if (fscanf(file, "%lf", &A[i]) != 1) {
            fprintf(stderr, "Failed to read A at index %d\n", i);
            free(A);
            free(b);
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    if (fscanf(file, " b =") != 0) {
        fprintf(stderr, "Failed to read 'b ='\n");
        free(A);
        free(b);
        fclose(file);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < params.N; i++) {
        if (fscanf(file, "%lf", &b[i]) != 1) {
            fprintf(stderr, "Failed to read b at index %d\n", i);
            free(A);
            free(b);
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
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

void generate_problem(double *A, double *b, CGParams params) {
    srand(SEED);
    if (A == NULL) {
        fprintf(stderr, "Memory allocation failed for A\n");
        exit(EXIT_FAILURE);
    }

    if (b == NULL) {
        fprintf(stderr, "Memory allocation failed for b\n");
        exit(EXIT_FAILURE);
    }
    generateRandomSPDMatrix(A, params.N);
    
    for (int i = 0; i < params.N; i++) {
        b[i] = (double)(rand() % 100);
    }
}

void print_cg_params(const CGParams params) {
    printf("Conjugate Gradient Parameters:\n");
    printf("Dimension N: %d\n", params.N);
    printf("Tolerance (tol): %g\n", params.tol);
    printf("Maximum Iterations (max_iter): %d\n", params.max_iter);
    printf("\n");
}