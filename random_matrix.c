#include "random_matrix.h"
#include <stdlib.h>
#include <time.h>

void generateRandomSparseMatrix(double matrix[N][N], int size, float sparsity) {
    int i, j;
    int nonZeroCount = (int)(size * size * (1 - sparsity)); // Total number of non-zero elements

    // Seed the random number generator
    srand((unsigned)time(NULL));

    // Initialize the matrix with zeros
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            matrix[i][j] = 0.0;
        }
    }

    // Randomly fill the matrix with non-zero values
    while (nonZeroCount > 0) {
        int row = rand() % size;
        int col = rand() % size;

        // Only add a non-zero element if the position is currently zero
        if (matrix[row][col] == 0) {
            matrix[row][col] = (double)rand() / RAND_MAX * 100.0; // Random value between 0 and 100
            nonZeroCount--;
        }
    }
}


// Function to generate a 1D vector
void generateVector(int vector[], int size) {
    for (int i = 0; i < size; i++) {
        vector[i] = rand() % 10; // Random integer between 0 and 9
    }
}


int countNonZeroElements(double matrix[N][N], int size) {
    int count = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (matrix[i][j] != 0.0) {
                count++;
            }
        }
    }
    return count;
}

 void printMatrix(const double matrix[N][N], int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%5.2f ", matrix[i][j]);
        }
        printf("\n");
    }
}
  
