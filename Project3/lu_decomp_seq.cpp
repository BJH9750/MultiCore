#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#define ADDR(A, N, X, Y) &A[X*N + Y]

void decomp(double *U, double *L, int n);
void row_calculation(double *A_pivot, double *A_target, double *LU_pivot, int size);
void matrix_print(double *mat, int size);

int main(int argc, char** argv){
    int n, s;
    double *A, *LU;
    n = atoi(argv[1]);
    s = atoi(argv[2]);

    srand(s);
    A = (double *)malloc(sizeof(double) * n * n);
    LU = (double *)malloc(sizeof(double) * n * n);
    double *L = (double *)malloc(sizeof(double) * n * n);

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i*n + j] = (double) (rand() % 100);
            LU[i*n + j] = A[i*n + j];
        }
    }

    decomp(LU, L, n);
    matrix_print(A, n);
    matrix_print(LU, n);
    matrix_print(L, n);
}

void decomp(double *U, double *L, int n){
    int pivot;
    for(pivot = 0; pivot < n - 1; ++pivot){
        double *U_pivot = ADDR(U, n, pivot, pivot);
        L[pivot*n + pivot] = 1;
        for(int next_p = pivot+1; next_p < n; ++next_p){
            double *U_target = ADDR(U, n, next_p, pivot);
            double *L_target = ADDR(L, n, next_p, pivot);
            row_calculation(U_pivot, U_target, L_target, n - pivot);
        }
    }
    L[pivot*n + pivot] = 1;
}

void row_calculation(double *A_pivot, double *A_target, double *L_target, int size){
    if(A_target[0] == 0) return;
    double coeff = A_target[0] / A_pivot[0];
    A_target[0] = 0;
    L_target[0] = coeff;
    for(int i = 1; i < size; ++i){
        A_target[i] = A_target[i] - coeff * A_pivot[i];
    }
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