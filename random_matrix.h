#ifndef RANDOM_MATRIX_H
#define RANDOM_MATRIX_H

#define N  10// Size of the matrix (N x N)
#define SPARSITY 0.90 // % of the matrix will be zeros

void generateRandomSparseMatrix(double matrix[N][N], int size, float sparsity);
int countNonZeroElements(double matrix[N][N], int size);
void printMatrix(const double matrix[N][N], int size); // Declare the print function
void generateVector(int vector[N], int size); // New function declaration

#endif // RANDOM_MATRIX_H
