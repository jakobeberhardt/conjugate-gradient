#ifndef MISC_H
#define MISC_H

typedef struct CGParams {
    double* A;   
    double* b;    
    double tol;   
    int max_iter; 
    int N;        
} CGParams;

void print_result(double *x, double *r, double *p, double *Ap, double residual, int iter, int N);
CGParams init_cg(const char* filename);
void print_cg_params(const CGParams params);
CGParams random_init_cg(const char* filename);
void generateRandomSPDMatrix(double *A, int N);
double wtime();

#endif /* MISC_H */