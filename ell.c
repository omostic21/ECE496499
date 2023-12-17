#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "random_matrix.h"
// Define structure for the ELL format
typedef struct Sparse_ELL {
    size_t nrows;
    size_t ncols;
    size_t max_nnz;
    size_t* indices;
    double** data;
} Sparse_ELL;

// Function to create a sparse matrix in the ELL format
int create_sparse_ell(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t max_nnz,
    Sparse_ELL* A_ell
);

// Function to print the ELL format representation
void print_sparse_ell(const Sparse_ELL* A_ell);

// Function to perform matrix-vector multiplication for an ELL format matrix
void matrix_vector_sparse_ell(const Sparse_ELL* A_ell, const int* vec, double* result);

int main() {
   double matrix[N][N]; // Use this 2D array to store the sparse matrix
    generateRandomSparseMatrix(matrix, N, SPARSITY);
    printMatrix(matrix, N);
    size_t n_rows = N;
    size_t n_cols = N;
    size_t max_nnz = countNonZeroElements(matrix, N);; // Set this to the maximum number of non-zero elements in any row

 int x[N];
    generateVector(x, N); // This will fill 'vector' with random values
    printf("\n\n");
    
for (int cn=0; cn<N; cn++){
    printf("%d\n", x[cn]);
}

 double Ax[N]; // Resultant vector


    // Create a sparse matrix in the ELL format
    Sparse_ELL A_ell;
    create_sparse_ell(matrix, n_rows, n_cols, max_nnz, &A_ell);

    // Print the ELL format representation
    print_sparse_ell(&A_ell);
    
    //Next 3 lines to declare timing metrics to measure CPU time
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    // Call multiplication function
    matrix_vector_sparse_ell(&A_ell, x, Ax);

     end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nResult: \n");
    for (size_t i = 0; i < n_rows; ++i) {
        printf("%2.2f\n", Ax[i]);
    }
    
    printf("\n");
    printf("Execution time: %f seconds\n", cpu_time_used);

    // Free memory
    for (size_t i = 0; i < A_ell.nrows; ++i) {
        free(A_ell.data[i]);
    }
    free(A_ell.data);
    free(A_ell.indices);

    return EXIT_SUCCESS;
}

int create_sparse_ell(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t max_nnz,
    Sparse_ELL* A_ell
) {
    A_ell->nrows = nrows;
    A_ell->ncols = ncols;
    A_ell->max_nnz = max_nnz;

    A_ell->indices = (size_t*)calloc(nrows * max_nnz, sizeof(size_t));
    A_ell->data = (double**)malloc(nrows * sizeof(double*));

    for (size_t i = 0; i < nrows; ++i) {
        A_ell->data[i] = (double*)calloc(max_nnz, sizeof(double));
    }

    for (size_t i = 0; i < nrows; ++i) {
        size_t nnz = 0;
        for (size_t j = 0; j < ncols; ++j) {
            if (A[i * ncols + j] != 0.0) {
                A_ell->data[i][nnz] = A[i * ncols + j];
                A_ell->indices[i * max_nnz + nnz] = j;
                nnz++;
            }
        }
    }

    return EXIT_SUCCESS;
}

void print_sparse_ell(const Sparse_ELL* A_ell) {
    printf("Number of rows: %zu\n", A_ell->nrows);
    printf("Number of columns: %zu\n", A_ell->ncols);
    printf("Maximum non-zzeros per row: %zu\n", A_ell->max_nnz);
    printf("Indices and Data:\n");

    for (size_t i = 0; i < A_ell->nrows; ++i) {
        printf("Row %zu: [", i);
        for (size_t j = 0; j < A_ell->max_nnz; ++j) {
            printf("(%zu, %2.2f) ", A_ell->indices[i * A_ell->max_nnz + j], A_ell->data[i][j]);
        }
        printf("]\n");
    }
}

void matrix_vector_sparse_ell(const Sparse_ELL* A_ell, const int* vec, double* result) {
    for (size_t i = 0; i < A_ell->nrows; ++i) {
        result[i] = 0.0;
        for (size_t j = 0; j < A_ell->max_nnz; ++j) {
            size_t col = A_ell->indices[i * A_ell->max_nnz + j];
            double value = A_ell->data[i][j];
            result[i] += value * vec[col];
        }
    }
}
