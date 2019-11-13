#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <algorithm>
#include <functional>

using namespace std;

#define ROOT 0
#define ADDR(A, N, X, Y) &A[X*N + Y]

void decomp(double *LU, int n);
void row_calculation(double *A_pivot, double *A_target, double *LU_pivot, int size);
double verify(double *A , double *LU);
void matrix_print(double *mat, int size);

int main(int argc, char** argv){
    int n, s;
    double *A, *LU;
    n = atoi(argv[1]);
    s = atoi(argv[2]);

    srand(s);
    A = (double *)malloc(sizeof(double) * n * n);
    LU = (double *)malloc(sizeof(double) * n * n);

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i*n + j] = (double) (rand() % 100);
            LU[i*n + j] = A[i*n + j];
        }
    }
    
    decomp(LU, n);
    matrix_print(A, n);
    matrix_print(LU, n);
    verify(A, LU);
}

void decomp(double *LU, int n){

    int size, rank;
    for(int pivot = 0; pivot < n - 1; ++pivot){
        double *LU_pivot = ADDR(LU, n, pivot, pivot);
        for(int next_p = pivot+1; next_p < n; ++next_p){
            double *LU_target = ADDR(LU, n, next_p, pivot);
            row_calculation(LU_pivot, LU_target, LU_target, n - pivot);
            //matrix_print(LU, n);
        }
    }
}

void row_calculation(double *A_pivot, double *A_target, double *LU_pivot, int size){
    double coeff = A_target[0] / A_pivot[0];
    LU_pivot[0] = coeff;
    for(int i = 1; i < size; ++i){
        LU_pivot[i] = A_target[i] - coeff * A_pivot[i];
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