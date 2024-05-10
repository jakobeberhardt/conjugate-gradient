#include <math.h>

int conjugate_gradient(const double *A, double *Ap, const double *b, double *x, double *r, double *p, double *rsold, const double *tol, double *residual, int nrow, int ii, int iter) {
	
	double sum, alpha, beta, rsnew, p_Ap, local_residual;

	// init R
	#pragma omp parallel for private(sum)
	for(int i=0; i<nrow; i++) {
		sum=0;
		for(int j=0; j<nrow; j++) {
			sum += A[i * nrow + j] * x[j];
		}
		r[i] = b[i] - sum;
	}

	// init p
	if(iter == 0) {
		for(int i=0; i<nrow; i++)
			p[i] = r[i];


		*rsold = 0.0;

		for(int i=0; i<nrow; i++)
			*rsold += r[i] * r[i];
	}

	// Compute Ap
	#pragma omp parallel for private(sum)
	for(int i=0; i<nrow; i++) {
		sum = 0;
		for(int j=0; j<nrow; j++) {
			sum += A[i * nrow + j] * p[j];
		}
		Ap[i] = sum;
	}

	// calculate dot product p Ap
	p_Ap = 0;
	for(int i=0; i<nrow; i++)
		p_Ap += p[i] * Ap[i];

	// Calculate alpha
	alpha = *rsold / p_Ap; 

	// Set x to xk + alpha*p and r to r -alpha*ap
	for(int i=0; i<nrow; i++) {
		x[i] = x[i] + (alpha * p[i]);
		r[i] = r[i] - (alpha*Ap[i]);
	}

	local_residual = 0;
	for(int i=0; i<nrow; i++)
    	local_residual += r[i] * r[i];

	*residual = sqrt(local_residual);
	if (*residual < *tol) {
    	return 1; 
	}

	// compute rk+1T * rk+1 as rsnew
	rsnew = 0.0;
	for(int i=0; i<nrow; i++)
		rsnew += r[i] * r[i];
	
	// compute beta rsnew / rsold
	beta = rsnew / *rsold;
	*rsold = rsnew;

	// update direction p
	for(int i=0; i<nrow; i++)
		p[i] = r[i] + (beta * p[i]);

	return 0;
}