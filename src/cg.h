#ifndef CG_H
#define CG_H

#include <stddef.h>

int conjugate_gradient(const double *A, const double *b, double *Ap, double *x, double *r, double *p, double *rsold, const double *tol, double *residual, int nrow, const int local_n, int const iter, int rank);

#endif /* CG_H */
