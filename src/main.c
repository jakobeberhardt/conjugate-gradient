#include <stdlib.h>
#include "cg.h"
#include "misc.h"
#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char **argv) {
	if (argc < 2) {
        fprintf(stderr, "Usage: %s [-r] <config file>\n", argv[0]);
        return 1;
    }

	int rank, size, iter,conv,randomInit, local_nrow, local_elements;
	double *A, *b, *x,*r, *p, *Ap, *local_A, *local_b, residual, rsold, start_time, end_time, flops_per_iteration, total_flops, flops_rate, gflops_rate;

	MPI_Init( &argc, &argv );
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	
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
		params = read_cg_config(filename);
		print_cg_params(params);
	} 
	
	MPI_Bcast(&params.N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&params.max_iter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&params.tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	local_nrow = getRowCount(params.N, size, rank);
	local_elements = (params.N * params.N ) / size;
	
	if(rank == 0){
		if (randomInit) {
		    A = (double*) malloc(params.N * params.N * sizeof(double));
			b = (double*) malloc(params.N * sizeof(double));
			x = calloc(params.N, sizeof(double));
			r = (double*)aligned_alloc(32, params.N * sizeof(double));
			p = calloc(params.N, sizeof(double));
			Ap = calloc(params.N, sizeof(double));
			local_A = (double*) malloc(local_elements * sizeof(double));
			local_b = (double*) malloc(local_nrow * sizeof(double));
			generate_problem(A, b, params);
		} else {
			load_problem(filename, A, b, params);
		}
	}
	
	if(rank !=0) {
		local_A = (double*) malloc(local_elements * sizeof(double));
		local_b = (double*) malloc(local_nrow * sizeof(double));
		x = calloc(params.N, sizeof(double));
		r = (double*)aligned_alloc(32, local_nrow * sizeof(double));
		p = calloc(params.N, sizeof(double));
		Ap = calloc(local_nrow, sizeof(double));
	}
	
	MPI_Scatter(A, local_elements, MPI_DOUBLE, local_A, local_elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(b, local_nrow, MPI_DOUBLE, local_b, local_nrow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	iter = 0;
	residual = 0.0;
	rsold = 0;
	conv = 0;
	if(rank == 0) start_time = wtime();
	do {
	 	conv = conjugate_gradient(local_A, local_b, Ap, x, r, p, &rsold, &params.tol, &residual, params.N, local_nrow, iter, rank);
		iter++;
		if(rank == 0)printf("Residual after iteration %d: %.17f\n", iter, residual);
	} while (iter <= params.max_iter && !conv);

	if(rank == 0) end_time = wtime();
	

	MPI_Gather((rank == 0 ? MPI_IN_PLACE : x), local_nrow, MPI_DOUBLE, x, local_nrow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather((rank == 0 ? MPI_IN_PLACE : p), local_nrow, MPI_DOUBLE, p, local_nrow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather((rank == 0 ? MPI_IN_PLACE : r), local_nrow, MPI_DOUBLE, r, local_nrow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gather((rank == 0 ? MPI_IN_PLACE : Ap), local_nrow, MPI_DOUBLE, Ap, local_nrow, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank == 0) {
		flops_per_iteration = 4 * params.N * params.N + 14 * params.N - 8;
		total_flops = flops_per_iteration * iter;
		flops_rate = total_flops / (end_time - start_time);
		gflops_rate = flops_rate / 1e9;
		
		printf("\n");
		print_result(x, r, p, Ap, residual,iter, params.N, conv);
		printf("\n");
		printf("Elapsed time   : %f seconds\n", end_time - start_time);
		printf("Iteration      : %d\n", iter);
		printf("Total FLOPs    : %e\n", total_flops);
		printf("GFLOPs Rate    : %.3f GFLOP/s\n", gflops_rate);
	}

	if(rank == 0) {
		free(A);
		free(b);
	}
	free(local_A);
	free(local_b);
	free(x);
	free(r);
	free(p);
	free(Ap);
	MPI_Finalize();
    return 0;
}