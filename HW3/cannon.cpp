#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <mpi.h>
#include <cmath>

#define P
int calculateBlock(int world, int n){
    world = sqrt(world);
    while((n % world) != 0){
        --world;
    }
    return world;
}

void printMatrix(double *A, char* str, int n, int rank){
    printf("-- rank : %d %s --\n", rank, str);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            printf("%3d ", (int)(A[i*n + j]));
        }
        printf("\n");
    }
    printf("------------\n");
}

void matrixProduct(double *A, double *B, double *C, int n){
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            for(int k = 0; k < n; ++k){
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }
}

int main(int argc, char** argv){
    int n = atoi(argv[1]);
    int seed = atoi(argv[2]);
    int rank, size;
    srand(seed);
    
    double *A, *B, *C;
    
    A = (double *)malloc(sizeof(double) * n * n);
    B = (double *)malloc(sizeof(double) * n * n);
    C = (double *)calloc(sizeof(double), n * n);

    int nblock, nproc, all_proc;
    //printf("%d\n", calculate_block(80, 600));

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    nproc = calculateBlock(size, n);
    nblock = n / nproc; // row, col num per block
    all_proc = nproc * nproc;
    if(rank == 0){
        for(int i = 0; i < n; ++i){
            for(int j = 0; j < n; ++j){
                A[i*n+j] = rand() % 100;
                B[i*n+j] = rand() % 100;
            }
        }
    }
    
    int color = (rank < all_proc) ? 0 : 1;
    MPI_Comm square_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &square_comm);

    if(color == 0){
        int nprocdim[] = {nproc, nproc};
        int period[] = {1, 1};
        int coords[2];
        int right=0, left=0, down=0, up=0;
        //printf("%d before comm rank\n", rank);
        MPI_Comm cart_comm;
        MPI_Cart_create(square_comm, 2, nprocdim, period, 1, &cart_comm);
        MPI_Comm_rank(cart_comm, &rank);
        //printf("%d comm rank\n", rank);

        int start[2] = {0, 0};
        int subsizes[2]  = {nblock, nblock};
        int bigsizes[2]  = {n, n};
        MPI_Datatype array_t, sub_array_t;
        MPI_Type_create_subarray(2, bigsizes, subsizes, start, MPI_ORDER_C, MPI_DOUBLE, &array_t);
        MPI_Type_create_resized(array_t, 0, nblock * sizeof(double), &sub_array_t);
        MPI_Type_commit(&sub_array_t);

        int *sendcounts = (int *)malloc(sizeof(int) * n);
        int *displs = (int *)malloc(sizeof(int) * n);

        if(rank == 0){
            for(int i = 0; i < n; ++i){
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

        if(rank == 0){
            printf("nproc : %d\tnblock : %d\n", nproc, nblock);
            for(int i = 0; i < nproc; ++i){
                for(int j = 0; j < nproc; ++j){
                    printf("%d ", displs[i*nproc + j]);
                }
            }
            printf("\n");
        }

        double *subA = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subB = (double *)malloc(sizeof(double) * nblock * nblock);
        double *subC = (double *)calloc(nblock * nblock, sizeof(double));

        MPI_Scatterv(A, sendcounts, displs, sub_array_t, subA, nblock * nblock, MPI_DOUBLE, 0, square_comm);
        MPI_Scatterv(B, sendcounts, displs, sub_array_t, subB, nblock * nblock, MPI_DOUBLE, 0, square_comm);
        //printf("%d scatter\n", rank);
        MPI_Cart_coords(cart_comm, rank, 2, coords);
        MPI_Cart_shift(cart_comm, 1, coords[0], &left, &right);
        MPI_Sendrecv_replace(subA, nblock * nblock, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
        MPI_Cart_shift(cart_comm, 0, coords[1], &up, &down);
        MPI_Sendrecv_replace(subB, nblock * nblock, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);

        for(int i = 0; i < nproc; ++i){
            matrixProduct(subA, subB, subC, nblock);

            MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
            MPI_Sendrecv_replace(subA, nblock * nblock, MPI_DOUBLE, left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);
            MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
            MPI_Sendrecv_replace(subB, nblock * nblock, MPI_DOUBLE, up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
        }

        MPI_Gatherv(subC, nblock * nblock, MPI_DOUBLE, C, sendcounts, displs, sub_array_t, 0, square_comm);

        #ifdef P
        if(rank == 0){
            printMatrix(A, "A", n, rank);
            printMatrix(B, "B", n, rank);
            printMatrix(C, "C", n, rank);
        }
        #endif
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}