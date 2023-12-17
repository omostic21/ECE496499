#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "random_matrix.h"

//Define structure for Coordiante Sparse row format
typedef struct Sparse_CSR {
    size_t n_rows;
    size_t n_cols;
    size_t n_nz;
    size_t* row_ptrs;
    size_t* col_indices;
    double* values;
} Sparse_CSR;

//Pointer memory assignment
int create_sparse_csr(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    size_t n_nz,
    Sparse_CSR* A_csr
);

int print_sparse_csr(const Sparse_CSR* A_csr);

int matrix_vector_sparse_csr(
    const Sparse_CSR* A_coo,
    const int* vec,
    double* res
);

int free_sparse_csr(Sparse_CSR* A_csr);


int main (int argc, char** argv) {
    double matrix[N][N]; // Use this 2D array to store the sparse matrix
    generateRandomSparseMatrix(matrix, N, SPARSITY);
    printMatrix(matrix, N);
    size_t n_rows = N;
    size_t n_cols = N;
    size_t n_nz = countNonZeroElements(matrix, N);

  
   printf("\n\n");

    int x[N];
    generateVector(x, N); // This will fill 'vector' with random values
    printf("\n\n");
    
for (int cn=0; cn<N; cn++){
    printf("%d\n", x[cn]);
}

    double Ax[N]; // Resultant vector

    // Start here- begin calling functions
    Sparse_CSR A_csr;

    create_sparse_csr(&matrix[0][0], n_rows, n_cols, n_nz, &A_csr);

    print_sparse_csr(&A_csr);

    //Next 3 lines to declare timing metrics to measure CPU time
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();
    
    matrix_vector_sparse_csr(&A_csr, x, Ax);

    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nResult: \n");
    
    for (size_t i=0; i<n_rows; ++i) {
        printf("%02.2f\n", Ax[i]);
    }

    printf("\nExecution time: %f seconds\n", cpu_time_used);  

    free_sparse_csr(&A_csr);

    return EXIT_SUCCESS;
}


int create_sparse_csr(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    size_t n_nz,
    Sparse_CSR* A_csr
) {
    A_csr->n_rows = n_rows;
    A_csr->n_cols = n_cols;
    A_csr->n_nz = n_nz;
    A_csr->row_ptrs = calloc(n_rows+1, sizeof(size_t));
    A_csr->col_indices = calloc(n_nz, sizeof(size_t));
    A_csr->values = calloc(n_nz, sizeof(double));

    size_t nz_id = 0;

    for (size_t i=0; i<n_rows; ++i) {
        A_csr->row_ptrs[i] = nz_id;
        for (size_t j=0; j<n_cols; ++j) {
            if (A[i*n_cols + j] != 0.0) {
                A_csr->col_indices[nz_id] = j;
                A_csr->values[nz_id] = A[i*n_cols + j];
                nz_id++;
            }
        }
    }

    A_csr->row_ptrs[n_rows] = nz_id;

    return EXIT_SUCCESS;
}

int print_sparse_csr(const Sparse_CSR* A_csr) {
    printf("row\tcol\tval\n");
    printf("----\n");
    for (size_t i=0; i<A_csr->n_rows; ++i) {
        size_t nz_start = A_csr->row_ptrs[i];
        size_t nz_end = A_csr->row_ptrs[i+1];
        for (size_t nz_id=nz_start; nz_id<nz_end; ++nz_id) {
            size_t j = A_csr->col_indices[nz_id];
            double val = A_csr->values[nz_id];
            printf("%d\t%d\t%02.2f\n", i, j, val);
        }
    }
    return EXIT_SUCCESS;
}

int matrix_vector_sparse_csr(
    const Sparse_CSR* A_csr,
    const int* vec,
    double* res
) {
    for (size_t i=0; i<A_csr->n_rows; ++i) {
        res[i] = 0.0;
        size_t nz_start = A_csr->row_ptrs[i];
        size_t nz_end = A_csr->row_ptrs[i+1];
        for (size_t nz_id=nz_start; nz_id<nz_end; ++nz_id) {
            size_t j = A_csr->col_indices[nz_id];
            double val = A_csr->values[nz_id];
            res[i] = res[i] + val * vec[j];
        }
    }
    return EXIT_SUCCESS;
}

int free_sparse_csr(Sparse_CSR* A_csr) {
    free(A_csr->row_ptrs);
    free(A_csr->col_indices);
    free(A_csr->values);

    return EXIT_SUCCESS;
}