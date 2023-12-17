#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 7 // Size of the A
#define SPARSITY 0.9 // 90% of the A will be zeros

void generateRandomSparseA(int A[N][N], int size, float sparsity) {
    int i, j;
    int nonZeroCount = (int)(size * size * (1 - sparsity));

    // Initialize the A with zeros
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            A[i][j] = 0;
        }
    }

    // Randomly fill the A with non-zero values
    while (nonZeroCount > 0) {
        int row = rand() % size;
        int col = rand() % size;

        if (A[row][col] == 0) {
            A[row][col] = rand() % 100 + 1; // Random value between 1 and 100
            nonZeroCount--;
        }
    }
}

int main() {
    int A[N][N];
    srand(time(NULL)); // Seed for random number generation

    generateRandomSparseA(A, N, SPARSITY);

    // Print the A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%d ", A[i][j]);
        }
        printf("\n");
    }

    return 0;
}
gcc -c random_matrix.c -o random_matrix.o
gcc -o csr_program csr.c random_matrix.o
./csr_program

gcc -c random_matrix.c -o random_matrix.o
gcc -o coo_program coo.c random_matrix.o
./coo_program

gcc -c random_matrix.c -o random_matrix.o
gcc -o ell_program ell.c random_matrix.o
./ell_program

gcc -c random_matrix.c -o random_matrix.o
gcc -o nocompprogram nocomp.c random_matrix.o
./nocompprogram