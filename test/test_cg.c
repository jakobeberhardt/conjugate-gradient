#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <math.h>
#include "mpi/mpi.h"
#include "cg.h"
#include "misc.h"

#define N 4
#define TOL 1e-6

void test_conjugate_gradient(void) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double A[N][N] = {
        {4.0, 1.0, 1.0, 0.0},
        {1.0, 3.0, 1.0, 0.0},
        {1.0, 1.0, 2.0, 0.0},
        {0.0, 0.0, 0.0, 1.0}
    };
    double b[] = {1.0, 2.0, 3.0, 4.0};
    double x[N] = {0.0};
    double r[N], p[N], Ap[N];
    double rsold = 0.0;
    double residual = 0.0;
    int max_iter = 1000;
    int iter = 0;
    int conv = 0;
    double tol = TOL;

    do {
        conv = conjugate_gradient(&A[0][0], Ap, b, x, r, p, &rsold, &tol, &residual, N, 0, iter, rank, size);
        iter++;
        MPI_Barrier(MPI_COMM_WORLD);
    } while (iter <= max_iter && !conv && rank == 0);

    if (rank == 0) {
        CU_ASSERT_DOUBLE_EQUAL(x[0], -0.176471, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(x[1], 0.235294, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(x[2], 1.470588, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(x[3], 4.000000, 1e-6);
        CU_ASSERT_DOUBLE_EQUAL(residual, 0.000000, 1e-6);
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    CU_initialize_registry();
    CU_pSuite suite = CU_add_suite("CG_Test_Suite", 0, 0);
    CU_add_test(suite, "test_conjugate_gradient", test_conjugate_gradient);
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    MPI_Finalize();
    return CU_get_error();
}
