#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#define INF 1e9

// Define structure for Coordiante Sparse row format
typedef struct Sparse_CSR {
    size_t nrows;
    size_t ncols;
    size_t n_nz;
    size_t* row_ptrs;
    size_t* col_indices;
    double* values;
} Sparse_CSR;

// Pointer memory assignment
int create_sparse_csr(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t n_nz,
    Sparse_CSR* A_csr
);

int print_sparse_csr(const Sparse_CSR* A_csr);

int matrix_vector_sparse_csr(
    const Sparse_CSR* A_coo,
    const double* vec,
    double* res
);

int free_sparse_csr(Sparse_CSR* A_csr);

void floyd_warshall_csr(const Sparse_CSR* A_csr, double** dist, size_t** next) {
    size_t n = A_csr->nrows;
    *dist = (double*)malloc(n * n * sizeof(double));
    *next = (size_t*)malloc(n * n * sizeof(size_t));

    // Initialize dist and next matrices
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                (*dist)[i * n + j] = 0;
            } else {
                (*dist)[i * n + j] = INF;
            }
            (*next)[i * n + j] = INF;
        }
    }

    // Fill dist matrix with the values from the CSR format
    for (size_t i = 0; i < n; ++i) {
        size_t nz_start = A_csr->row_ptrs[i];
        size_t nz_end = A_csr->row_ptrs[i + 1];
        for (size_t nz_id = nz_start; nz_id < nz_end; ++nz_id) {
            size_t j = A_csr->col_indices[nz_id];
            double val = A_csr->values[nz_id];
            (*dist)[i * n + j] = val;
            (*next)[i * n + j] = j;
        }
    }

    // Floyd-Warshall algorithm
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if ((*dist)[i * n + k] != INF && (*dist)[k * n + j] != INF &&
                    (*dist)[i * n + k] + (*dist)[k * n + j] < (*dist)[i * n + j]) {
                    (*dist)[i * n + j] = (*dist)[i * n + k] + (*dist)[k * n + j];
                    (*next)[i * n + j] = (*next)[i * n + k];
                }
            }
        }
    }
}

int main(int argc, char** argv) {
    size_t nrows = 5;
    size_t ncols = 5;
    size_t n_nz = 12;

    double A[] = {
        9, 0, 0, 16, 0,
        5, 4, 7, 5, 0,
        5, 0, 0, 8, 12,
        0, 0, 19, 5, 0,
        0, 0, 0, 0, 14
    };
    double x[] = {
        1,
        2,
        3,
        4,
        5
    };
    double Ax[5];

    // Start here - begin calling functions
    Sparse_CSR A_csr;

    create_sparse_csr(A, nrows, ncols, n_nz, &A_csr);

    print_sparse_csr(&A_csr);

    matrix_vector_sparse_csr(&A_csr, x, Ax);

    for (size_t i = 0; i < nrows; ++i) {
        printf("%02.2f\n", Ax[i]);
    }

    double* dist;
    size_t* next;
    floyd_warshall_csr(&A_csr, &dist, &next);

    printf("All-Pairs Shortest Paths (Distance Matrix):\n");
    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j) {
            if (dist[i * nrows + j] == INF) {
                printf("INF ");
            } else {
                printf("%2.2f ", dist[i * nrows + j]);
            }
        }
        printf("\n");
    }

    free_sparse_csr(&A_csr);
    free(dist);
    free(next);

    

    return EXIT_SUCCESS;
}

int create_sparse_csr(
    const double* A,
    size_t nrows,
    size_t ncols,
    size_t n_nz,
    Sparse_CSR* A_csr
) {
    A_csr->nrows = nrows;
    A_csr->ncols = ncols;
    A_csr->n_nz = n_nz;
    A_csr->row_ptrs = calloc(nrows + 1, sizeof(size_t));
    A_csr->col_indices = calloc(n_nz, sizeof(size_t));
    A_csr->values = calloc(n_nz, sizeof(double));

    size_t nz_id = 0;

    for (size_t i = 0; i < nrows; ++i) {
        A_csr->row_ptrs[i] = nz_id;
        for (size_t j = 0; j < ncols; ++j) {
            if (A[i * ncols + j] != 0.0) {
                A_csr->col_indices[nz_id] = j;
                A_csr->values[nz_id] = A[i * ncols + j];
                nz_id++;
            }
        }
    }

    A_csr->row_ptrs[nrows] = nz_id;

    return EXIT_SUCCESS;
}

int print_sparse_csr(const Sparse_CSR* A_csr) {
    printf("row\tcol\tval\n");
    printf("----\n");
    for (size_t i = 0; i < A_csr->nrows; ++i) {
        size_t nz_start = A_csr->row_ptrs[i];
        size_t nz_end = A_csr->row_ptrs[i + 1];
        for (size_t nz_id = nz_start; nz_id < nz_end; ++nz_id) {
            size_t j = A_csr->col_indices[nz_id];
            double val = A_csr->values[nz_id];
            printf("%zu\t%zu\t%02.2f\n", i, j, val);
        }
    }
    return EXIT_SUCCESS;
}

int matrix_vector_sparse_csr(
    const Sparse_CSR* A_csr,
    const double* vec,
    double* res
) {
    for (size_t i = 0; i < A_csr->nrows; ++i) {
        res[i] = 0.0;
        size_t nz_start = A_csr->row_ptrs[i];
        size_t nz_end = A_csr->row_ptrs[i + 1];
        for (size_t nz_id = nz_start; nz_id < nz_end; ++nz_id) {
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


