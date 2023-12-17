#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "random_matrix.h" // Make sure this is correctly implemented as discussed earlier

typedef struct Sparse_Hybrid {
   size_t nrows;
    size_t ncols;
    size_t max_nz_per_row; // Maximum number of non-zeros per row for ELL
    double** ell_data;
    size_t** ell_indices;
    size_t coo_nnz; // Number of non-zero elements for COO
    size_t* coo_row_indices;
    size_t* coo_col_indices;
    double* coo_values;
} Sparse_Hybrid;

int create_sparse_hybrid(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t max_nz_per_row,
    Sparse_Hybrid* A_hybrid
);

int print_sparse_hybrid(const Sparse_Hybrid* A_hybrid);

int matrix_vector_sparse_hybrid(
    const Sparse_Hybrid* A_hybrid,
    const int* vec,
    double* res
);

int main() {
    size_t max_nz_per_row = 3; // Adjust based on maximum expected non-zeros per row
    double matrix[N][N];
    generateRandomSparseMatrix(matrix, N, SPARSITY);
    printMatrix(matrix, N); // Print the generated random matrix

    size_t n_rows = N;
    size_t n_cols = N;
    size_t n_nz = countNonZeroElements(matrix, N); // Count non-zero elements

    int x[N]; // Integer vector
    generateVector(x, N); // Fill 'x' with random integer values
    for (int cn = 0; cn < N; cn++) {
        printf("%d ", x[cn]);
    }
    printf("\n");

    double Ax[N]; // Resultant vector

    Sparse_Hybrid A_hybrid;
    create_sparse_hybrid(&matrix[0][0], n_rows, n_cols, max_nz_per_row, &A_hybrid);

   // Timing metrics
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();


    matrix_vector_sparse_hybrid(&A_hybrid, x, Ax); // Perform matrix-vector multiplication

    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    
    printf("\nResult:\n");
    for (size_t i = 0; i < n_rows; ++i) {
        printf("%2.2f\n", Ax[i]);
    }
    free_sparse_hybrid(&A_hybrid); // Free the memory allocated for the hybrid matrix

    return 0;
}


int create_sparse_hybrid(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t max_nz_per_row,
    Sparse_Hybrid* A_hybrid
) {
    A_hybrid->nrows = nrows;
    A_hybrid->ncols = ncols;
    A_hybrid->max_nz_per_row = max_nz_per_row;
    A_hybrid->ell_data = (double**)malloc(nrows * sizeof(double*));
    A_hybrid->ell_indices = (size_t**)malloc(nrows * sizeof(size_t*));
    A_hybrid->coo_nnz = 0; // Initialize COO non-zero count

    // Allocate ELL data
    for (size_t i = 0; i < nrows; ++i) {
        A_hybrid->ell_data[i] = (double*)calloc(max_nz_per_row, sizeof(double));
        A_hybrid->ell_indices[i] = (size_t*)calloc(max_nz_per_row, sizeof(size_t));
    }

    // Populate ELL and COO portions
    size_t* nnz_per_row = calloc(nrows, sizeof(size_t)); // Track non-zeros per row
    for (size_t row = 0; row < nrows; ++row) {
        for (size_t col = 0; col < ncols; ++col) {
            double val = A[row * ncols + col];
            if (val != 0.0) {
                // Check if we can put it in ELL part
                if (nnz_per_row[row] < max_nz_per_row) {
                    A_hybrid->ell_data[row][nnz_per_row[row]] = val;
                    A_hybrid->ell_indices[row][nnz_per_row[row]] = col;
                    nnz_per_row[row]++;
                } else {
                    // Put it in the COO part
                    A_hybrid->coo_row_indices = (size_t*)realloc(A_hybrid->coo_row_indices, (A_hybrid->coo_nnz + 1) * sizeof(size_t));
                    A_hybrid->coo_col_indices = (size_t*)realloc(A_hybrid->coo_col_indices, (A_hybrid->coo_nnz + 1) * sizeof(size_t));
                    A_hybrid->coo_values = (double*)realloc(A_hybrid->coo_values, (A_hybrid->coo_nnz + 1) * sizeof(double));
                    A_hybrid->coo_row_indices[A_hybrid->coo_nnz] = row;
                    A_hybrid->coo_col_indices[A_hybrid->coo_nnz] = col;
                    A_hybrid->coo_values[A_hybrid->coo_nnz] = val;
                    A_hybrid->coo_nnz++;
                }
            }
        }
    }
    free(nnz_per_row); // Free temporary array

    return 0;
}

int print_sparse_hybrid(const Sparse_Hybrid* A_hybrid) {
    printf("Hybrid format\n");
    printf("nrows: %zu, ncols: %zu, max_nz_per_row: %zu, coo_nnz: %zu\n", 
           A_hybrid->nrows, A_hybrid->ncols, A_hybrid->max_nz_per_row, A_hybrid->coo_nnz);

    printf("ELL format:\n");
    for (size_t i = 0; i < A_hybrid->nrows; ++i) {
        printf("Row %zu: ", i);
        for (size_t j = 0; j < A_hybrid->max_nz_per_row; ++j) {
            printf("(%zu, %2.2f) ", A_hybrid->ell_indices[i][j], A_hybrid->ell_data[i][j]);
        }
        printf("\n");
    }

    printf("COO format:\n");
    for (size_t i = 0; i < A_hybrid->coo_nnz; ++i) {
        printf("(%zu, %zu, %2.2f) ", A_hybrid->coo_row_indices[i], A_hybrid->coo_col_indices[i], A_hybrid->coo_values[i]);
    }
    printf("\n");

    return 0;
}

int matrix_vector_sparse_hybrid(
    const Sparse_Hybrid* A_hybrid,
    const int* vec,
    double* res
) {
    // Initialize result vector
    for (size_t i = 0; i < A_hybrid->nrows; ++i) {
        res[i] = 0.0;
    }

    // Multiply ELL part
    for (size_t i = 0; i < A_hybrid->nrows; ++i) {
        for (size_t j = 0; j < A_hybrid->max_nz_per_row; ++j) {
            if (A_hybrid->ell_indices[i][j] < A_hybrid->ncols) { // Valid column index
                size_t col = A_hybrid->ell_indices[i][j];
                double val = A_hybrid->ell_data[i][j];
                res[i] += val * vec[col];
            }
        }
    }

    // Multiply COO part
    for (size_t i = 0; i < A_hybrid->coo_nnz; ++i) {
        size_t row = A_hybrid->coo_row_indices[i];
        size_t col = A_hybrid->coo_col_indices[i];
        double val = A_hybrid->coo_values[i];
        res[row] += val * vec[col];
    }

    return 0;
}

int free_sparse_hybrid(Sparse_Hybrid* A_hybrid) {
    for (size_t i = 0; i < A_hybrid->nrows; ++i) {
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
