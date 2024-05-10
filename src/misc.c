#include <stdio.h>

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