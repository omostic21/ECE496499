#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "random_matrix.h"

typedef struct coo {
    size_t n_rows;
    size_t n_cols;
    size_t nnz;
    size_t* row_indices;
    size_t* col_indices;
    double* values;
} coo;


int create_coo(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    size_t nnz,
    coo* A_coo
);

int print_coo(const coo* A_coo);

int matrix_vector_coo(
    const coo* A_coo,
    const int* vec,
    double* res
);

int free_coo(coo* A_coo);


int main (int argc, char** argv) {
    double matrix[N][N]; // Use this 2D array to store the sparse matrix
    generateRandomSparseMatrix(matrix, N, SPARSITY);
   printMatrix(matrix, N);
    size_t n_rows = N;
    size_t n_cols = N;
    size_t nnz = countNonZeroElements(matrix, N);

  
    printf("\n\n");

    int x[N];
    generateVector(x, N); // This will fill 'vector' with random values
    printf("\n\n");
    
for (int cn=0; cn<N; cn++){
    printf("%d\n", x[cn]);
}

    double Ax[N]; // Resultant vector

   coo A_coo;

   create_coo(matrix, n_rows, n_cols, nnz, &A_coo);

   print_coo(&A_coo);

   //Next 3 lines to declare timing metrics to measure CPU time
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    matrix_vector_coo(&A_coo, x, Ax);
    
    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("\nResult:\n");
    for (size_t i=0; i<n_cols; ++i) {
        printf("%2.2f\n", Ax[i]);
    }

 printf("\nExecution time: %f seconds\n", cpu_time_used);  

    free_coo(&A_coo);

    return EXIT_SUCCESS;
}


int create_coo(
    const double* A,
    size_t n_rows,
    size_t n_cols,
    size_t nnz,
    coo* A_coo
) {
    A_coo->n_rows = n_rows;
    A_coo->n_cols = n_cols;
    A_coo->nnz = nnz;
    A_coo->row_indices = calloc(nnz, sizeof(size_t));
    A_coo->col_indices = calloc(nnz, sizeof(size_t));
    A_coo->values = calloc(nnz, sizeof(double));

    size_t nnz_id = 0;

    for (size_t i=0; i<n_rows; ++i) {
        for (size_t j=0; j<n_cols; ++j) {
            if (A[i*n_cols + j] != 0) {
                A_coo->row_indices[nnz_id] = i;
                A_coo->col_indices[nnz_id] = j;
                A_coo->values[nnz_id] = A[i*n_cols + j];
                nnz_id++;
            }
        }
    }

    return EXIT_SUCCESS;
}

int print_coo(const coo* A_coo) {
    printf("\n");
    printf("row\tcol\tval\n");
    printf("---\t---\t---\n");
    for(size_t nnz_id=0; nnz_id<A_coo->nnz; ++nnz_id) {
        size_t row_id = A_coo->row_indices[nnz_id];
        size_t col_id = A_coo->col_indices[nnz_id];
        double value = A_coo->values[nnz_id];

        printf("%d\t%d\t%02.2f\n", row_id, col_id, value);
    }

    return EXIT_SUCCESS;
}

int matrix_vector_coo(
    const coo* A_coo,
    const int* vec,
    double* res
) {
    for (size_t i=0; i<A_coo->n_cols; ++i) {
        res[i] = 0.0;
    }

    for (size_t nnz_id=0; nnz_id<A_coo->nnz; ++nnz_id) {
        size_t row_id = A_coo->row_indices[nnz_id];
        size_t col_id = A_coo->col_indices[nnz_id];
        double value = A_coo->values[nnz_id];

        res[row_id] += value * vec[col_id];
    }

    return EXIT_SUCCESS;
}

int free_coo(coo* A_coo) {
    free(A_coo->row_indices);
    free(A_coo->col_indices);
    free(A_coo->values);

    return EXIT_SUCCESS;
}
