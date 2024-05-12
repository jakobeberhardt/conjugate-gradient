#ifndef CG_H
#define CG_H

#include <stddef.h>

int conjugate_gradient(const double *A, double *Ap, const double *b, double *x, double *r, double *p, double *rsold, const double *tol, double *residual, int nrow, const int ii, int const iter, int rank, int size);

#endif /* CG_H */
