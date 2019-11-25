#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>

#define ROOT 0
#define ADDR(A, N, X, Y) &A[X*N + Y]

using namespace std;

void decomp(double *A, double *U, double *L, int n, int rank, int size);
void subDecomp(double *U, int n, int rank, int size);
void rowCalculation(double *A_pivot, double *A_target, double *L_target, int size);
void split_matrix(double *U, double *L, int n);
double* matrixProduct(double *A, double *U, double *L, int size, int n, int rank);
void matrix_print(double *mat, int size);
int calculateBlock(int world, int n);
void subProduct(double *A, double *B, double *C, int n, int reset);
double verify(double *A, double *LU, int size, int n, int rank);
double *matrixIeverse(double *A, int n);
void subSubtract(double *A, double *B, double *C, int n);
void cleanMatrix(double *A, int n, char c);


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
    double *A, *tmpA, *U, *L, *LU = NULL;
    int rank, size;
    n = atoi(argv[1]);
    s = atoi(argv[2]);

    srand(s);
    A = (double *)malloc(sizeof(double) * n * n);
    tmpA = (double *)malloc(sizeof(double) * n * n);
    U = (double *)malloc(sizeof(double) * n * n);
    L = (double *)malloc(sizeof(double) * n * n);
    double starttime, endtime;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            A[i*n + j] = (int) (rand() % 1000);
            tmpA[i*n + j] = A[i*n + j];
        }
    }

    decomp(tmpA, U, L, n, rank, size);
    if(rank == 0) printf("decomp\n");
    cleanMatrix(U, n, 'U');
    cleanMatrix(L, n, 'L');
    if(rank == 0) printf("clear\n");
    LU = matrixProduct(A, U, L, size, n, rank);
    if(rank == 0) printf("product\n");
    double value = verify(A, LU, size, n, rank); 

    if(rank == 0){
        cout << "verify : "  << setprecision(15) << value << endl;
    } 

    MPI_Finalize();
}

void setIdentity(double *A, int n){
    memset(A, 0, sizeof(double)* n * n);
    for(int i = 0; i < n; ++i){
        A[i*n + i] = 1;
    }
}

double* matrixIeverse(double *A, int n){
    double *B = (double *)malloc(sizeof(double) * n * n);
    setIdentity(B, n);
    double *tmpA = (double *)malloc(sizeof(double) * n * n);
    memcpy(tmpA, A, sizeof(double) * n * n);
    int row;
    double pivot;
    for(int i = 0; i < n; ++i){
        row = i * n;
        pivot = tmpA[row + i];
        for(int j = 0; j < n; ++j){
            tmpA[row + j] /= pivot;
            B[row + j] /= pivot;
        }
        int x;
        for(int j = 0, x = (i+1)%n; j < n-1; ++j, x = (x+1)%n){
            pivot = tmpA[x*n + i];
            for(int k = 0; k < n; ++k){
                if(k >= i) tmpA[x*n + k] -= pivot * tmpA[row + k];
                B[x*n + k] -= pivot * B[row + k];
            }
        }
    }
    return B;
}

void subDecomp(double *U, double *L, int n){
    int pivot;
    for(pivot = 0; pivot < n - 1; ++pivot){
        double *U_pivot = ADDR(U, n, pivot, pivot);
        L[pivot*n + pivot] = 1;
        for(int next_p = pivot+1; next_p < n; ++next_p){
            double *U_target = ADDR(U, n, next_p, pivot);
            double *L_target = ADDR(L, n, next_p, pivot);
            rowCalculation(U_pivot, U_target, L_target, n - pivot);
        }
    }
    L[pivot*n + pivot] = 1;
}

void rowCalculation(double *A_pivot, double *A_target, double *L_target, int size){
    if(A_target[0] == 0) return;
    double coeff = A_target[0] / A_pivot[0];
    A_target[0] = 0;
    L_target[0] = coeff;
    for(int i = 1; i < size; ++i){
        A_target[i] = A_target[i] - coeff * A_pivot[i];
    }
}

void decomp(double *A, double *U, double *L, int n, int rank, int size){
    
    int nblock, nproc, all_proc;
    nproc = calculateBlock(size, n);
    nblock = n / nproc; // row, col num per block
    all_proc = nproc * nproc;
    int color = (rank < all_proc) ? 0 : 1;
    MPI_Comm square_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &square_comm);

    if(color == 0){
        int nprocdim[] = {nproc, nproc};
        int period[] = {1, 1};
        int coords[2];
        int right=0, left=0, down=0, up=0;

        MPI_Comm cart_comm, sub_cart_comm;
        MPI_Cart_create(square_comm, 2, nprocdim, period, 1, &cart_comm);
        MPI_Comm_rank(cart_comm, &rank);
        MPI_Cart_coords(cart_comm, rank, 2, coords);

        int start[2] = {0, 0};
        int subsizes[2]  = {nblock, nblock};
        int bigsizes[2]  = {n, n};
        MPI_Datatype array_t, sub_array_t;
        MPI_Type_create_subarray(2, bigsizes, subsizes, start, MPI_ORDER_C, MPI_DOUBLE, &array_t);
        MPI_Type_create_resized(array_t, 0, nblock * sizeof(double), &sub_array_t);
        MPI_Type_commit(&sub_array_t);

        int *sendcounts = (int *)malloc(sizeof(int) * all_proc);
        int *displs = (int *)malloc(sizeof(int) * all_proc);

        MPI_Comm sub_comm, root_comm;
        int sub_rank;
        int sub_size;
        int sub_coords[2];

        double *subA = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subL = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subU = (double *)malloc(sizeof(double) * nblock * nblock);
        double *invL, *invU;

        int diag_rank;
        int diag_color = (coords[0] == coords[1]) ? 0 : 1;
        MPI_Comm_split(square_comm, diag_color, rank, &root_comm);
        MPI_Comm_rank(root_comm, &diag_rank);
        
        for(int m = 0; m < nproc; ++m){
            int sub_color = (coords[0] >= m && coords[1] >= m) ? 0 : 1;
            MPI_Comm_split(square_comm, sub_color, rank, &sub_comm);
            if(sub_color == 0){
                MPI_Cart_create(sub_comm, 2, nprocdim, period, 1, &sub_cart_comm);
                MPI_Comm_rank(sub_comm, &sub_rank);
                MPI_Comm_size(sub_comm, &sub_size);
                MPI_Cart_coords(sub_cart_comm, rank, 2, sub_coords);

                sendcounts = (int *)malloc(sizeof(int) * sub_size);
                displs = (int *)malloc(sizeof(int) * sub_size);

                nprocdim[0] -= 1;
                nprocdim[1] -= 1;

                if(sub_rank == 0){
                    for(int i = 0; i < sub_size; ++i){
                        sendcounts[i] = 1;
                    }
                    int displ = 0;
                    for(int i = 0; i < nproc - m; ++i){
                        for(int j = 0; j < nproc - m; ++j, ++displ){
                            displs[i*(nproc - m) + j] = displ; 
                        }
                        displ += (nblock - 1) * nproc + (nblock - 1) * m;
                    }
                }
                if(rank == 6) printf("asd\n");
                MPI_Scatterv(ADDR(A, n, m*nblock, m*nblock), sendcounts, displs, sub_array_t, subA, nblock * nblock, MPI_DOUBLE, ROOT, sub_comm);

                if(sub_rank == ROOT){
                    subDecomp(subA, subL, nblock);
                    memcpy(subU, subA, sizeof(double) * nblock * nblock);
                    cleanMatrix(subU, nblock, 'U');
                    cleanMatrix(subL, nblock, 'L');
                }

                if(sub_size != 1){
                    MPI_Comm row_comm, col_comm;
                    int row_color = sub_rank / (nproc - m);
                    int col_color = sub_rank % (nproc - m);

                    MPI_Comm_split(sub_comm, row_color, sub_rank, &row_comm);
                    MPI_Comm_split(sub_comm, col_color, sub_rank, &col_comm);
                    
                    int row_rank, row_size, col_rank, col_size;
                    MPI_Comm_rank(row_comm, &row_rank);
                    MPI_Comm_size(row_comm, &row_size);
                    MPI_Comm_rank(col_comm, &col_rank);
                    MPI_Comm_size(col_comm, &col_size);

                    if(row_color == 0){
                        MPI_Bcast(subL, nblock * nblock, MPI_DOUBLE, ROOT, row_comm);
                        if(sub_rank != ROOT){
                            invL = matrixIeverse(subL, nblock);
                            subProduct(invL, subA, subU, nblock, 1);
                        }
                    }
                    
                    if(col_color == 0){
                        MPI_Bcast(subU, nblock * nblock, MPI_DOUBLE, ROOT, col_comm);
                        if(sub_rank != ROOT){
                            invU = matrixIeverse(subU, nblock);
                            subProduct(subA, invU, subL, nblock, 1);
                        }
                    }

                    if(row_color != 0){
                        int row_root;
                        int row_root_coords[] = {row_color, 0};
                        MPI_Cart_rank(cart_comm, row_root_coords, &row_root);
                        MPI_Bcast(subL, nblock * nblock, MPI_DOUBLE, ROOT, row_comm);
                    }

                    if(col_color != 0){
                        int col_root;
                        int col_root_coords[] = {0, col_color};
                        MPI_Cart_rank(cart_comm, col_root_coords, &col_root);
                        MPI_Bcast(subU, nblock * nblock, MPI_DOUBLE, ROOT, col_comm);
                    }

                    if(row_color != 0 && col_color != 0){
                        double *tmpA = (double *)malloc(sizeof(double) * nblock * nblock);
                        subProduct(subL, subU, tmpA, nblock, 1);
                        subSubtract(subA, tmpA, subA, nblock);
                    }
                }

                MPI_Gatherv(subA, nblock * nblock, MPI_DOUBLE, ADDR(A, n, m*nblock, m*nblock), sendcounts, displs, sub_array_t, ROOT, sub_comm);
                MPI_Gatherv(subL, nblock * nblock, MPI_DOUBLE, ADDR(L, n, m*nblock, m*nblock), sendcounts, displs, sub_array_t, ROOT, sub_comm);
                MPI_Gatherv(subU, nblock * nblock, MPI_DOUBLE, ADDR(U, n, m*nblock, m*nblock), sendcounts, displs, sub_array_t, ROOT, sub_comm);
                
                MPI_Comm_free(&sub_cart_comm);
                
                free(sendcounts);
                free(displs);
            }
            
            if(diag_color == 0){
                MPI_Bcast(A, n * n, MPI_DOUBLE, m, root_comm);
                MPI_Bcast(L, n * n, MPI_DOUBLE, m, root_comm);
                MPI_Bcast(U, n * n, MPI_DOUBLE, m, root_comm);
            }
            MPI_Barrier(square_comm); 
            MPI_Comm_free(&sub_comm);
        }
        MPI_Type_free(&sub_array_t);
        MPI_Type_free(&array_t);
        MPI_Comm_free(&root_comm);
        MPI_Comm_free(&cart_comm);
    }
    MPI_Comm_free(&square_comm);
    MPI_Barrier(MPI_COMM_WORLD);    

}

void cleanMatrix(double *A, int n, char c){
    int row, i, j;
    if(c == 'U'){
        for(i = 1; i < n; ++i){
            row = i * n;
            for(j = 0; j < i; ++j){
                A[row + j] = 0;
            }
        }
    }else{
        for(i = 0; i < n; ++i){
            row = i * n;
            for(j = i+1; j < n; ++j){
                A[row + j] = 0;
            }
        }
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
    MPI_Comm product_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &product_comm);
    
    if(color == 0){
        LU = (double *)calloc(sizeof(double), n * n);
        int nprocdim[] = {nproc, nproc};
        int period[] = {1, 1};
        int coords[2];
        int right=0, left=0, down=0, up=0;

        MPI_Comm cart_comm;
        MPI_Cart_create(product_comm, 2, nprocdim, period, 1, &cart_comm);
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

        double *subL = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subU = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subLU = (double *)calloc(nblock * nblock, sizeof(double));

        MPI_Scatterv(L, sendcounts, displs, sub_array_t, subL, nblock * nblock, MPI_DOUBLE, 0, product_comm);
        MPI_Scatterv(U, sendcounts, displs, sub_array_t, subU, nblock * nblock, MPI_DOUBLE, 0, product_comm);

        MPI_Cart_coords(cart_comm, rank, 2, coords);
        MPI_Cart_shift(cart_comm, 1, coords[0], &left, &right);
        MPI_Sendrecv_replace(subL, nblock * nblock, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
        MPI_Cart_shift(cart_comm, 0, coords[1], &up, &down);
        MPI_Sendrecv_replace(subU, nblock * nblock, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);

        for(int i = 0; i < nproc; ++i){
            subProduct(subL, subU, subLU, nblock, 0);

            MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
            MPI_Sendrecv_replace(subL, nblock * nblock, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
            MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
            MPI_Sendrecv_replace(subU, nblock * nblock, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
        }

        MPI_Gatherv(subLU, nblock * nblock, MPI_DOUBLE, LU, sendcounts, displs, sub_array_t, 0, product_comm);
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

void subProduct(double *A, double *B, double *C, int n, int reset){
    if(reset == 1) memset(C, 0, sizeof(double) * n * n);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

void subSubtract(double *A, double *B, double *C, int n){
    int row;
    for(int i = 0; i < n; ++i){
        row = i * n;
        for(int j = 0; j < n; ++j){
            C[row + j] = A[row + j] - B[row + j];
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
            row_sum = 0;
            for(int j = 0; j < n; ++j){
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
            row_sum = 0;
            if(size == 1){
                if(rank == 0){
                    for(int j = 0; j < n; ++j){

                        row_sum += (A[idx+j] - LU[idx+j]) * (A[idx+j] - LU[idx+j]);
                    }
                    row_sum = sqrt(row_sum);
                    sum += row_sum;
                }
            }else{
                int color = rank < size ? 0 : 1;
                MPI_Comm_split(MPI_COMM_WORLD, color, rank, &row_root);
                if(color == 0){
                    MPI_Scatter(A + idx, n, MPI_DOUBLE, A_row, n, MPI_DOUBLE, ROOT, row_root);
                    MPI_Scatter(LU + idx, n, MPI_DOUBLE, LU_row, n, MPI_DOUBLE, ROOT, row_root);

                    for(int j = 0; j < n; ++j){
                        row_sum += (A_row[j] - LU_row[j]) * (A_row[j] - LU_row[j]);
                    }

                    row_sum = sqrt(row_sum);

                    MPI_Reduce(&row_sum, &tmp_sum, 1, MPI_DOUBLE, MPI_SUM, ROOT, row_root);
                    sum += tmp_sum;
                }
            }
            
        }
    }

    return sum;
}
