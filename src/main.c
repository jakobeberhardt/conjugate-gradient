#include <stdlib.h>
#include "cg.h"
#include "misc.h"

#define N 4

int main(int argc, char **argv) {
double A[N][N] = {
        {4.0, 1.0, 1.0, 0.0},
        {1.0, 3.0, 1.0, 0.0},
        {1.0, 1.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}
    };
	double b[] = {1.0, 2.0, 3.0, 4.0};
	double *x = calloc(N, sizeof(double));
	double *r = calloc(N, sizeof(double));
	double *p = calloc(N, sizeof(double));
	double *Ap = calloc(N, sizeof(double));

	double tol = 1e-6;
	int max_iter = 1000;
	int iter = 0;
	double residual = 0.0;
	double rsold = 0;
	int conv = 0;

	do {
	 	conv = conjugate_gradient(&A[0][0], Ap, b, x, r, p, &rsold, &tol, &residual, N, 0, iter);
		iter++;

	} while (iter <= max_iter && !conv);

	print_result(x, r, p, Ap, residual,iter, N);

	free(x);
    free(r);
    free(p);
    free(Ap);

    return 0;
}