#include <stdlib.h>
#include "cg.h"
#include "misc.h"
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "mpi/mpi.h"

int main(int argc, char **argv) {
	if (argc < 2) {
        fprintf(stderr, "Usage: %s [-r] <config file>\n", argv[0]);
        return 1;
    }

    
	MPI_Init( &argc, &argv );
	int rank, size, iter,conv,randomInit;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	double *x,*r, *p, *Ap, residual, rsold, start_time, end_time;
	CGParams params;
    randomInit = 0;
    char *filename;

	if(rank == 0) {
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

		x = calloc(params.N, sizeof(double));
		r = calloc(params.N, sizeof(double));
		p = calloc(params.N, sizeof(double));
		Ap = calloc(params.N, sizeof(double));
	
	} 

	iter = 0;
	residual = 0.0;
	rsold = 0;
	conv = 0;

	MPI_Bcast(&params.N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&params.max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&params.tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	
	if(rank != 0) {
		params.A = (double*) malloc((params.N * params.N * sizeof(double)) / size);
		params.b = (double*) malloc((params.N * sizeof(double)) / size);
		x = calloc(params.N / size, sizeof(double));
		r = calloc(params.N / size, sizeof(double));
		p = calloc(params.N / size, sizeof(double));
		Ap = calloc(params.N / size, sizeof(double));
		printf("Allocated\n");
	}
	
	if(rank == 0) start_time = wtime();
	do {
	 	conv = conjugate_gradient(params.A, Ap, params.b, x, r, p, &rsold, &params.tol, &residual, params.N, 0, iter, rank, size);
		if(rank == 0)printf("Residual in iteration %d: %.17f\n", iter, residual);
		iter++;

	} while (iter <= params.max_iter && !conv);


	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0) {

		end_time = wtime();

		double flops_per_iteration = 4 * params.N * params.N + 14 * params.N - 8;
		double total_flops = flops_per_iteration * iter;
		double flops_rate = total_flops / (end_time - start_time);
		double gflops_rate = flops_rate / 1e9;
		
		printf("\n");
		print_result(x, r, p, Ap, residual,iter, params.N);
		printf("\n");
		printf("Elapsed time   : %f seconds\n", end_time - start_time);
		printf("Iteration      : %d\n", iter);
		printf("Total FLOPs    : %e\n", total_flops);
		printf("GFLOPs Rate    : %.3f GFLOP/s\n", gflops_rate);
	}
	free(params.A);
	free(params.b);
	free(x);
	free(r);
	free(p);
	free(Ap);
	MPI_Finalize();
    return 0;
}