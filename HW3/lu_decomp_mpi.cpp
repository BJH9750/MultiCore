#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;

#define ROOT 0
#define ADDR(A, N, X, Y) &A[X*N + Y]

void decomp(double *U, int n, int rank, int size);
void row_calculation(double *A_pivot, double *A_target, int size, int rank);
double verify(double *A , double *LU);
void matrix_print(double *mat, int size);

void print_row(double *A, int size, int rank, int sender){
    if(rank !=sender) printf("rank %d sender %d", rank, sender);
    else printf("rank %d calculate", rank);
    for(int i = 0; i < size; ++i){
        printf("%6.3lf ", A[i]);
    }
    printf("\n");
}

int main(int argc, char** argv){
    int n, s;
    double *A, *U;
    int rank, size;
    n = atoi(argv[1]);
    s = atoi(argv[2]);

    srand(s);
    A = (double *)malloc(sizeof(double) * n * n);
    U = (double *)malloc(sizeof(double) * n * n);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i*n + j] = (double) (rand() % 100);
            U[i*n + j] = A[i*n + j];
        }
    }

    
    decomp(U, n, rank, size);

    if(rank == 0){
        matrix_print(A, n);
        matrix_print(U, n);
    }

    verify(A, U);
    MPI_Finalize();
}

void decomp(double *U, int n, int rank, int size){
    for(int pivot = 0; pivot < n - 1; ++pivot){
        double *U_pivot = ADDR(U, n, pivot, pivot);
        for(int next_p = pivot+1; next_p < n; ++next_p){
            if(next_p % size == rank){
                double *U_target = ADDR(U, n, next_p, pivot);
                row_calculation(U_pivot, U_target, n - pivot, rank);
            }
        }
        for(int next_p = pivot+1; next_p < n; ++next_p){
            double *U_target = ADDR(U, n, next_p, pivot);
            MPI_Bcast(U_target, n - pivot, MPI_DOUBLE, next_p % size, MPI_COMM_WORLD);
        }
    }
}

void row_calculation(double *A_pivot, double *A_target, int size, int rank){
    //print_row(A_target, size, rank, rank);
    if(A_target[0] == 0) return;
    double coeff = A_target[0] / A_pivot[0];
    A_target[0] = coeff;
    for(int i = 1; i < size; ++i){
        A_target[i] = A_target[i] - coeff * A_pivot[i];
    }
}

double verify(double *A, double *LU){
    return 0;
}

void matrix_print(double *mat, int size){
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            printf("%6.3lf\t", mat[i*size+j]);
        }
        printf("\n");
    }
    printf("\n");
}