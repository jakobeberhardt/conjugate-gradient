#include <math.h>
#include <mpi/mpi.h>
#include <misc.h>

int conjugate_gradient(const double *A, double *Ap, const double *b, double *x, double *r, double *p, double *rsold, const double *tol, double *residual, int N, int ii, int iter, int rank, int size) {
	double sum, alpha, beta, rsnew, local_rsnew, p_Ap, local_residual, local_rsold, local_p_Ap;
	int local_n = N / size;

	if(iter == 0) {
		// init R
		for(int i=0; i<local_n; i++) {
			sum=0;
			for(int j=0; j<N; j++) {
				sum += A[i * N + j] * x[j];
			}
			r[i] = b[i] - sum;
		}

		// init p
		for(int i=0; i<local_n; i++)
			p[i + rank * local_n] = r[i];

		local_rsold = 0.0;
		for(int i=0; i<local_n; i++)
			local_rsold += r[i] * r[i];

		MPI_Allreduce(&local_rsold, rsold, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	}

	// share direction p
    MPI_Allgather(p + rank * local_n, local_n, MPI_DOUBLE, p, local_n, MPI_DOUBLE, MPI_COMM_WORLD);

	// Compute Ap
	for(int i=0; i<local_n; i++) {
		sum = 0.0;
		for(int j=0; j<N; j++) {
			sum += A[i * N + j] * p[j];
		}
		Ap[i] = sum;
	}

	// calculate dot product p Ap
	p_Ap = 0.0;
	local_p_Ap = 0.0;
	for(int i=0; i<local_n; i++) {
		local_p_Ap += p[i + rank * local_n] * Ap[i];
	}

	MPI_Allreduce(&local_p_Ap, &p_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	// Calculate alpha
	alpha = *rsold / p_Ap; 

	// Set x to xk + alpha*p and r to r -alpha*ap
	for(int i=0; i<local_n; i++) {
		x[i] = x[i] + (alpha * p[i + rank * local_n]);
		r[i] = r[i] - (alpha*Ap[i]);
	}

	// rk+1T * rk+1 as rsnew
	rsnew = 0.0;
	local_rsnew = 0.0;
	for(int i=0; i<local_n; i++)
    	local_rsnew += r[i] * r[i];

	MPI_Allreduce(&local_rsnew, &rsnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	*residual = sqrt(rsnew);
	if (*residual < *tol) {
    	return 1; 
	}

	// compute beta rsnew / rsold
	beta = rsnew / *rsold;
	*rsold = rsnew;

	// update direction p
	for(int i=0; i<local_n; i++)
		p[i + rank * local_n] = r[i] + (beta * p[i + rank * local_n]);
		
	return 0;
}