#ifndef MISC_H
#define MISC_H

typedef struct CGParams {
    double tol;   
    int max_iter; 
    int N;        
} CGParams;

void print_result(double *x, double *r, double *p, double *Ap, double residual, int iter, int N, int conv);
CGParams read_cg_config(const char* filename);
void print_cg_params(const CGParams params);
void load_problem(const char* filename, double *A, double *b, CGParams params);
void generate_problem(double *A, double *b, CGParams params);
double wtime();
int getRowCount(const int totalRows, const int size, const int rank);

#endif /* MISC_H */