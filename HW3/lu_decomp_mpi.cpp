#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#define ROOT 0
#define ADDR(A, N, X, Y) &A[X*N + Y]

using namespace std;

void decomp(double *U, int n, int rank, int size);
void row_calculation(double *A_pivot, double *A_target, int size, int rank);
void split_matrix(double *U, double *L, int n);
double* matrixProduct(double *A, double *U, double *L, int size, int n, int rank);
void matrix_print(double *mat, int size);
int calculateBlock(int world, int n);
void subProduct(double *A, double *B, double *C, int n);
double verify(double *A, double *LU, int size, int n, int rank);

MPI_Comm square_comm;

void print_row(double *A, int size, int rank, int sender){
    if(rank !=sender) printf("rank %d sender %d", rank, sender);
    else printf("rank %d calculate", rank);
    for(int i = 0; i < size; ++i){
        printf("%lf ", A[i]);
    }
    printf("\n");
}

int main(int argc, char** argv){
    int n, s;
    double *A, *U, *L, *LU = NULL;
    int rank, size;
    n = atoi(argv[1]);
    s = atoi(argv[2]);

    srand(s);
    A = (double *)malloc(sizeof(double) * n * n);
    U = (double *)malloc(sizeof(double) * n * n);
    L = (double *)malloc(sizeof(double) * n * n);

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i*n + j] = (int) (rand() % 1000);
            U[i*n + j] = A[i*n + j];
            if(i == j) L[i*n + j] = 1;
        }
    }

    decomp(U, n, rank, size);
    if(rank == 0 && argc > 3) matrix_print(A, n);
    split_matrix(U, L, n);
    LU = matrixProduct(A, U, L, size, n, rank);
    if(rank == 0 && argc > 3) matrix_print(LU, n);
    double value = verify(A, LU, size, n, rank); 
    if(rank == 0){
        cout << "verify : "  << setprecision(15) << value << endl;
    } 

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

void split_matrix(double *U, double *L, int n){
    for(int i = 1; i < n; ++i){
        memcpy(ADDR(L, n, i, 0), ADDR(U, n, i, 0), sizeof(double) * i);
        memset(ADDR(U, n, i, 0), 0, sizeof(double) * i);
    }
}

double* matrixProduct(double *A, double *U, double *L, int size, int n, int rank){

    int nblock, nproc, all_proc;
    double *LU;
    nproc = calculateBlock(size, n);
    nblock = n / nproc; // row, col num per block
    all_proc = nproc * nproc;

    int color = (rank < all_proc) ? 0 : 1;

    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &square_comm);

    if(color == 0){
        LU = (double *)calloc(sizeof(double), n * n);
        int nprocdim[] = {nproc, nproc};
        int period[] = {1, 1};
        int coords[2];
        int right=0, left=0, down=0, up=0;

        MPI_Comm cart_comm;
        MPI_Cart_create(square_comm, 2, nprocdim, period, 1, &cart_comm);
        MPI_Comm_rank(cart_comm, &rank);

        int start[2] = {0, 0};
        int subsizes[2]  = {nblock, nblock};
        int bigsizes[2]  = {n, n};
        MPI_Datatype array_t, sub_array_t;
        MPI_Type_create_subarray(2, bigsizes, subsizes, start, MPI_ORDER_C, MPI_DOUBLE, &array_t);
        MPI_Type_create_resized(array_t, 0, nblock * sizeof(double), &sub_array_t);
        MPI_Type_commit(&sub_array_t);

        int *sendcounts = (int *)malloc(sizeof(int) * all_proc);
        int *displs = (int *)malloc(sizeof(int) * all_proc);

        if(rank == 0){
            for(int i = 0; i < all_proc; ++i){
                sendcounts[i] = 1;
            }
            int displ = 0;
            for(int i = 0; i < nproc; ++i){
                for(int j = 0; j < nproc; ++j, ++displ){
                    displs[i*nproc + j] = displ; 
                }
                displ += (nblock-1) * nproc;
            }
        }
        #ifdef PA
        if(rank == 0){
            printf("nproc : %d\tnblock : %d\n", nproc, nblock);
            for(int i = 0; i < nproc; ++i){
                for(int j = 0; j < nproc; ++j){
                    printf("%d ", displs[i*nproc + j]);
                }
            }
            printf("\n");
        }
        #endif
        double *subL = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subU = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subLU = (double *)calloc(nblock * nblock, sizeof(double));

        MPI_Scatterv(L, sendcounts, displs, sub_array_t, subL, nblock * nblock, MPI_DOUBLE, 0, square_comm);
        MPI_Scatterv(U, sendcounts, displs, sub_array_t, subU, nblock * nblock, MPI_DOUBLE, 0, square_comm);
        //printf("%d scatter\n", rank);
        MPI_Cart_coords(cart_comm, rank, 2, coords);
        MPI_Cart_shift(cart_comm, 1, coords[0], &left, &right);
        MPI_Sendrecv_replace(subL, nblock * nblock, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
        MPI_Cart_shift(cart_comm, 0, coords[1], &up, &down);
        MPI_Sendrecv_replace(subU, nblock * nblock, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);

        for(int i = 0; i < nproc; ++i){
            subProduct(subL, subU, subLU, nblock);

            MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
            MPI_Sendrecv_replace(subL, nblock * nblock, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
            MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
            MPI_Sendrecv_replace(subU, nblock * nblock, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
        }

        MPI_Gatherv(subLU, nblock * nblock, MPI_DOUBLE, LU, sendcounts, displs, sub_array_t, 0, square_comm);
        #ifdef PB
        if(rank == ROOT){
            matrix_print(U, n);
            matrix_print(L, n);
            matrix_print(A, n);
            matrix_print(LU, n);
        }
        #endif
        MPI_Type_free(&sub_array_t);
    }

    MPI_Barrier(MPI_COMM_WORLD);    

    return LU;
}

void matrix_print(double *mat, int size){
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            printf("%lf\t", mat[i*size+j]);
        }
        printf("\n");
    }
    printf("\n");
}


int calculateBlock(int world, int n){
    world = sqrt(world);
    world = world > n ? sqrt(n) : world;
    while((n % world) != 0){
        --world;
    }
    return world;
}

void subProduct(double *A, double *B, double *C, int n){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

double verify(double *A, double *LU, int size, int n, int rank){
    double *A_row = (double *)malloc(sizeof(double) * n);
    double *LU_row = (double *)malloc(sizeof(double) * n);

    double row_sum = 0;
    double sum = 0;
    double tmp_sum=0;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm row_root;

    if(size >= n){
        int color = rank < n ? 0 : 1;
        MPI_Comm_split(MPI_COMM_WORLD, color, rank, &row_root);
        if(color == 0){
            MPI_Scatter(A, n, MPI_DOUBLE, A_row, n, MPI_DOUBLE, ROOT, row_root);
            MPI_Scatter(LU, n, MPI_DOUBLE, LU_row, n, MPI_DOUBLE, ROOT, row_root);
            //printf("rank :%d size: %d\n", rank, size);
            for(int i = 0 ; i < n; ++i){
                //printf("%6.lf ", A_row[i]);
                row_sum += (A_row[i] - LU_row[i]) * (A_row[i] - LU_row[i]);
            }
            //printf("\n");
            row_sum = sqrt(row_sum);

            MPI_Reduce(&row_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, ROOT, row_root);
        }
    }else{
        int loop = n / size;
        int i;

        for(i = 0; i < loop; ++i){
            MPI_Scatter(A + i * n * size, n, MPI_DOUBLE, A_row, n, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            MPI_Scatter(LU + i * n * size, n, MPI_DOUBLE, LU_row, n, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
            #ifdef P
                printf("color 0 i : %d rank :%d size: %d\n", i, rank, size);
            #endif
            //if(rank == ROOT) printf("i : %d loop : %d size : %d\n", i, loop, size);
            row_sum = 0;
            for(int j = 0; j < n; ++j){
                #ifdef P
                    printf("%f ", A_row[j]);
                #endif
                row_sum += (A_row[j] - LU_row[j]) * (A_row[j] - LU_row[j]);
            }
            #ifdef P
                printf("\n");
            #endif
            row_sum = sqrt(row_sum);

            MPI_Reduce(&row_sum, &tmp_sum, 1, MPI_DOUBLE, MPI_SUM, ROOT, MPI_COMM_WORLD);
            sum += tmp_sum;
        }
        if(loop * size != n){
            int idx = i * size * n;
            size = n - loop * size;
            //if(rank == ROOT) printf("loop out i : %d loop : %d size : %d\n", i, loop, size);
            row_sum = 0;
            if(size == 1){
                if(rank == 0){
                    #ifdef P
                        printf("loop out color 0 rank :%d size: %d\n", rank, size);
                    #endif
                    for(int j = 0; j < n; ++j){
                        #ifdef P
                            printf("%f ", A[idx + j]);
                        #endif
                        row_sum += (A[idx+j] - LU[idx+j]) * (A[idx+j] - LU[idx+j]);
                    }
                    #ifdef P
                        printf("\n");
                    #endif
                    row_sum = sqrt(row_sum);
                    sum += row_sum;
                }
            }else{
                int color = rank < size ? 0 : 1;
                MPI_Comm_split(MPI_COMM_WORLD, color, rank, &row_root);
                if(color == 0){
                    MPI_Scatter(A + idx, n, MPI_DOUBLE, A_row, n, MPI_DOUBLE, ROOT, row_root);
                    MPI_Scatter(LU + idx, n, MPI_DOUBLE, LU_row, n, MPI_DOUBLE, ROOT, row_root);
                    #ifdef P
                        if(color == 0) printf("loop out color 0 rank :%d size: %d\n", rank, size);
                    #endif
                    for(int j = 0; j < n; ++j){
                        #ifdef P
                            if(color == 0) printf("%f ", A_row[j]);
                        #endif
                        row_sum += (A_row[j] - LU_row[j]) * (A_row[j] - LU_row[j]);
                    }
                    #ifdef P
                        if(color == 0) printf("\n");
                    #endif
                    row_sum = sqrt(row_sum);

                    MPI_Reduce(&row_sum, &tmp_sum, 1, MPI_DOUBLE, MPI_SUM, ROOT, row_root);
                    sum += tmp_sum;
                }
            }
            
        }
    }

    return sum;
}
