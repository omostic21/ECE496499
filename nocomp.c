#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "random_matrix.h"

int matrix_vector_sparse(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    const int* vec,
    double* res
);

int main (int argc, char** argv) {
    double matrix[N][N]; // Use this 2D array to store the sparse matrix
    generateRandomSparseMatrix(matrix, N, SPARSITY);
    size_t n_rows = N;
    size_t n_cols = N;

    int x[N];
    generateVector(x, N);

    double Ax[N]; // Resultant vector

    // Next 3 lines to declare timing metrics to measure CPU time
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();
    
    matrix_vector_sparse(&matrix[0][0], n_rows, n_cols, x, Ax);

    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nResult: \n");
    /*
    for (size_t i=0; i<n_rows; ++i) {
        printf("%02.2f\n", Ax[i]);
    }
    */

    printf("\nExecution time: %f seconds\n", cpu_time_used);

    return EXIT_SUCCESS;
}

int matrix_vector_sparse(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    const int* vec,
    double* res
) {
    for (size_t i = 0; i < n_rows; ++i) {
        res[i] = 0.0;
        for (size_t j = 0; j < n_cols; ++j) {
            res[i] += A[i * n_cols + j] * vec[j];
        }
    }
    return EXIT_SUCCESS;
}
