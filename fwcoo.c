#include <stdio.h>
#include <stdlib.h>
#include <time.h>

typedef struct Sparse_Hybrid {
    size_t nrows;
    size_t ncols;
    size_t ell_max_nnz;
    double** ell_data;
    size_t** ell_indices;
    size_t coo_nnz;
    size_t* coo_row_indices;
    size_t* coo_col_indices;
    double* coo_values;
} Sparse_Hybrid;

int create_sparse_hybrid(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t ell_max_nnz,
    Sparse_Hybrid* A_hybrid
);

int print_sparse_hybrid(const Sparse_Hybrid* A_hybrid);

int matrix_vector_sparse_hybrid(
    const Sparse_Hybrid* A_hybrid,
    const double* vec,
    double* res
);

int free_sparse_hybrid(Sparse_Hybrid* A_hybrid);

int main() {
    size_t nrows = 5;
    size_t ncols = 5;
    size_t ell_max_nnz = 3;  // Adjust this to match your matrix

    double A[] = {
        23, 0, 0, 5, 0,
        1, 99, 3, 2, 0,
        6, 0, 0, 8, 22,
        0, 0, 49, 5, 0,
        0, 0, 0, 0, 16
    };
    double x[] = {
        1,
        2,
        3,
        4,
        5
    };
    double Ax[5];

    Sparse_Hybrid A_hybrid;
    create_sparse_hybrid(A, nrows, ncols, ell_max_nnz, &A_hybrid);

    print_sparse_hybrid(&A_hybrid);

    //Next 3 lines to declare timing metrics to measure CPU time
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();

    matrix_vector_sparse_hybrid(&A_hybrid, x, Ax);

    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Execution time: %f seconds\n", cpu_time_used);

    printf("Matrix-Vector Multiplication:\n");
    for (size_t i = 0; i < nrows; ++i) {
        printf("%2.2f\n", Ax[i]);
    }
   
    printf("\nExecution time: %f seconds\n", cpu_time_used);
   
    free_sparse_hybrid(&A_hybrid);

    return 0;
}

int create_sparse_hybrid(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t ell_max_nnz,
    Sparse_Hybrid* A_hybrid
) {
    A_hybrid->nrows = nrows;
    A_hybrid->ncols = ncols;
    A_hybrid->ell_max_nnz = ell_max_nnz;
    A_hybrid->ell_data = (double**)malloc(ell_max_nnz * sizeof(double*));
    A_hybrid->ell_indices = (size_t**)malloc(ell_max_nnz * sizeof(size_t*));
    A_hybrid->coo_nnz = 0;

    for (size_t i = 0; i < ell_max_nnz; ++i) {
        A_hybrid->ell_data[i] = (double*)malloc(ncols * sizeof(double));
        A_hybrid->ell_indices[i] = (size_t*)malloc(ncols * sizeof(size_t));
    }

    for (size_t i = 0; i < ell_max_nnz; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            A_hybrid->ell_data[i][j] = 0.0;
            A_hybrid->ell_indices[i][j] = 0;
        }
    }

    for (size_t row = 0; row < nrows; ++row) {
        for (size_t col = 0; col < ncols; ++col) {
            double val = A[row * ncols + col];
            if (val != 0) {
                if (A_hybrid->coo_nnz == 0) {
                    A_hybrid->coo_row_indices = (size_t*)malloc(sizeof(size_t));
                    A_hybrid->coo_col_indices = (size_t*)malloc(sizeof(size_t));
                    A_hybrid->coo_values = (double*)malloc(sizeof(double));
                } else {
                    A_hybrid->coo_row_indices = (size_t*)realloc(A_hybrid->coo_row_indices, (A_hybrid->coo_nnz + 1) * sizeof(size_t));
                    A_hybrid->coo_col_indices = (size_t*)realloc(A_hybrid->coo_col_indices, (A_hybrid->coo_nnz + 1) * sizeof(size_t));
                    A_hybrid->coo_values = (double*)realloc(A_hybrid->coo_values, (A_hybrid->coo_nnz + 1) * sizeof(double));
                }
                A_hybrid->coo_row_indices[A_hybrid->coo_nnz] = row;
                A_hybrid->coo_col_indices[A_hybrid->coo_nnz] = col;
                A_hybrid->coo_values[A_hybrid->coo_nnz] = val;
                A_hybrid->coo_nnz++;
            }
        }
    }

    return 0;
}

int print_sparse_hybrid(const Sparse_Hybrid* A_hybrid) {
    printf("Hybrid format\n");
    printf("nrows: %zu, ncols: %zu, ell_max_nnz: %zu, coo_nnz: %zu\n", A_hybrid->nrows, A_hybrid->ncols, A_hybrid->ell_max_nnz, A_hybrid->coo_nnz);

    printf("ELL format:\n");
    for (size_t i = 0; i < A_hybrid->ell_max_nnz; ++i) {
        printf("Row %zu:\n", i);
        for (size_t j = 0; j < A_hybrid->ncols; ++j) {
            printf("  Data: %2.2f, Index: %zu\n", A_hybrid->ell_data[i][j], A_hybrid->ell_indices[i][j]);
        }
    }

    printf("COO format:\n");
    for (size_t i = 0; i < A_hybrid->coo_nnz; ++i) {
        printf("Row: %zu, Col: %zu, Value: %2.2f\n", A_hybrid->coo_row_indices[i], A_hybrid->coo_col_indices[i], A_hybrid->coo_values[i]);
    }

    return 0;
}

int matrix_vector_sparse_hybrid(
    const Sparse_Hybrid* A_hybrid,
    const double* vec,
    double* res
) {
    for (size_t i = 0; i < A_hybrid->nrows; ++i) {
        res[i] = 0.0;
    }

    for (size_t i = 0; i < A_hybrid->ell_max_nnz; ++i) {
        for (size_t j = 0; j < A_hybrid->ncols; ++j) {
            size_t row = A_hybrid->ell_indices[i][j];
            double val = A_hybrid->ell_data[i][j];
            res[row] += val * vec[j];
        }
    }

    for (size_t i = 0; i < A_hybrid->coo_nnz; ++i) {
        size_t row = A_hybrid->coo_row_indices[i];
        size_t col = A_hybrid->coo_col_indices[i];
        double val = A_hybrid->coo_values[i];
        res[row] += val * vec[col];
    }

    return 0;
}

int free_sparse_hybrid(Sparse_Hybrid* A_hybrid) {
    for (size_t i = 0; i < A_hybrid->ell_max_nnz; ++i) {
        free(A_hybrid->ell_data[i]);
        free(A_hybrid->ell_indices[i]);
    }
    free(A_hybrid->ell_data);
    free(A_hybrid->ell_indices);
    free(A_hybrid->coo_row_indices);
    free(A_hybrid->coo_col_indices);
    free(A_hybrid->coo_values);

    return 0;
}
