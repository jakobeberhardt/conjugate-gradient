#include <stdlib.h>
#include "cg.h"
#include "misc.h"
#include <stdio.h>
#include <time.h>
#include <string.h>

int main(int argc, char **argv) {
	if (argc < 2) {
        fprintf(stderr, "Usage: %s [-r] <config file>\n", argv[0]);
        return 1;
    }

    CGParams params;
    int randomInit = 0;
    char *filename;

    if (argc == 3 && strcmp(argv[1], "-r") == 0) {
        randomInit = 1;
        filename = argv[2];
    } else {
        filename = argv[1];
    }

    if (randomInit) {
        params = random_init_cg(filename);
    } else {
        params = init_cg(filename);
    }
	print_cg_params(params);

	double *x = calloc(params.N, sizeof(double));
	double *r = calloc(params.N, sizeof(double));
	double *p = calloc(params.N, sizeof(double));
	double *Ap = calloc(params.N, sizeof(double));
	int iter = 0;
	double residual = 0.0;
	double rsold = 0;
	int conv = 0;
	
	clock_t start, end;  
    double cpu_time_used;
	start = clock();

	do {
	 	conv = conjugate_gradient(params.A, Ap, params.b, x, r, p, &rsold, &params.tol, &residual, params.N, 0, iter);
		printf("Residual in iteration %d: %.17f\n", iter, residual);
		iter++;

	} while (iter <= params.max_iter && !conv);
	
	end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;  
    printf("Elapsed time   : %f seconds\n", cpu_time_used);
	printf("Iteration      : %d\n", iter);

	//print_result(x, r, p, Ap, residual,iter, params.N);

	free(x);
    free(r);
    free(p);
    free(Ap);

    return 0;
}