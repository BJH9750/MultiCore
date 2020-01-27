#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_BANKS 32
#define LOG_NUM_BANKS 5
#define CONFLICT_FREE_OFFSET(n) ((n) >> NUM_BANKS + (n) >> (2 * LOG_NUM_BANKS))

__global__ void countArray(int *array, int *histogram, int size, int max_val);
__global__ void scanCount(int *histogram, int *blockSum, int size);
__global__ void spreadSum(int *histogram, int *blockSum, int size);
__global__ void countSort(int *array, int *histogram, int size);
__host__ unsigned int maxTwo(int n);

__host__ void debug(const char* name, int *array, int *count, int *blockSum, int size, int max_val, int block){
    FILE* out = fopen(name, "w");
    if(array != NULL){
        for (int i = 0; i < size; ++i)
            fprintf(out, "array[%03d] = %d\n", i, array[i]);
        fprintf(out, "\n");
    }
    if(count != NULL){
        for (int i = 0; i < max_val; ++i)
            fprintf(out, "count[%03d] = %d\n", i, count[i]);
        fprintf(out, "\n");
    }
    if(blockSum != NULL){
        for (int i = 0; i < block; ++i){
            fprintf(out, "blockSum[%03d] = %d\n", i, blockSum[i]);
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

__host__ void printError(const char* loc){
    cudaError_t error = cudaGetLastError();
    printf("%s : %s\n", loc, cudaGetErrorString(error));
}

__host__ void counting_sort(int array[], int size, int max_val)
{
    int *d_array, *d_count, *d_blockSum;
    unsigned int nearMaxTwo = maxTwo(max_val);
    int block, thread = 1024;

    if(size <= 1024){
        thread = size;
        block = 1;
    }else{
        while((size % thread) != 0) --thread;
        block = size / thread;
    }

    int scanBlock, scanThread;
    scanThread = (nearMaxTwo > 1024) ? 1024 : nearMaxTwo;
    scanBlock = (nearMaxTwo == scanThread) ? 1 : (nearMaxTwo / scanThread) / 2;

    cudaMalloc(&d_array, sizeof(int) * size);
    cudaMalloc(&d_count, sizeof(int) * nearMaxTwo);
    cudaMalloc(&d_blockSum, sizeof(int) * scanBlock * 2);

    cudaMemcpy(d_array, array, sizeof(int) * size, cudaMemcpyHostToDevice);
    cudaMemset(d_count, 0, sizeof(int) * nearMaxTwo);
    cudaMemset(d_blockSum, 0, sizeof(int) * scanBlock * 2);
    
    countArray<<<block, thread>>>(d_array, d_count, size, max_val);

    scanCount<<<scanBlock, scanThread, scanThread * 2 * sizeof(int)>>>(d_count, d_blockSum, scanThread * 2);

    if(scanBlock <= 2048){
        scanCount<<<1, scanBlock / 2, scanBlock * sizeof(int)>>>(d_blockSum, NULL, scanBlock);
        spreadSum<<<scanBlock - 1, scanThread, 2 * scanThread * sizeof(int)>>>(d_count, d_blockSum, 2 * scanThread);

    }else{
        int scanBlock_r, scanThread_r, blockSumSize;
        int *d_blockSum_r[5] = {d_blockSum, NULL, NULL, NULL, NULL};
        int idx = 1;
        do{
            scanThread_r = 1024;
            scanBlock_r = (scanBlock / scanThread_r) / 2;
            blockSumSize = sizeof(int) * scanBlock_r * 2;
            cudaMalloc(&d_blockSum_r[idx], blockSumSize);
            cudaMemset(d_blockSum_r[idx], 0, blockSumSize);

            scanCount<<<scanBlock_r, scanThread_r, scanThread_r * 2 * sizeof(int)>>>(d_blockSum_r[idx - 1], d_blockSum_r[idx], scanThread_r * 2);

            ++idx;
        }while(scanBlock_r > 1024);
        --idx;

        scanCount<<<1, scanBlock_r / 2, scanBlock_r * sizeof(int)>>>(d_blockSum_r[idx], NULL, scanBlock_r);

        for(int i = idx; i > 0; --i){
            spreadSum<<<scanBlock_r - 1, scanThread_r, 2 * scanThread_r * sizeof(int)>>>(d_blockSum_r[idx-1], d_blockSum_r[idx], 2 * scanThread_r);
            scanBlock_r *= 2;
        }

        spreadSum<<<scanBlock - 1, scanThread, 2 * scanThread * sizeof(int)>>>(d_count, d_blockSum, 2 * scanThread);
        for(int i = 1; i < idx; ++i) cudaFree(d_blockSum_r[i]);
    }

    countSort<<<block, thread>>>(d_array, d_count, size);
    cudaMemcpy(array, d_array, sizeof(int) * size, cudaMemcpyDeviceToHost);

    cudaFree(d_array);
    cudaFree(d_count);
    cudaFree(d_blockSum);
}

__global__ void countArray(int *array, int *histogram, int size, int max_val){

    int taskPerBlock = (size > gridDim.x) ? (size / gridDim.x) : 1;
    int taskPerThread = (taskPerBlock > blockDim.x) ? (taskPerBlock / blockDim.x) : 1;

    int threadOffset = threadIdx.x + blockIdx.x * blockDim.x;
    int threadStart = threadOffset * taskPerThread;
    
    for(int i = threadStart; i < threadStart + taskPerThread; ++i){
        atomicAdd(&histogram[array[i]], 1);
    }
    __syncthreads();
}

__global__ void scanCount(int *histogram, int* blockSum, int size){
    extern __shared__ int buffer[];

    int threadId = threadIdx.x;
    int blockOffset = blockIdx.x * size;

    int offset = 1;
    int left = threadId;
    int right = threadId + (size / 2);
    int bankOffsetA = CONFLICT_FREE_OFFSET(left);
    int bankOffsetB = CONFLICT_FREE_OFFSET(left);
    buffer[left + bankOffsetA] = histogram[left + blockOffset];
    buffer[right + bankOffsetB] = histogram[right + blockOffset]; 
    int saveLeft = buffer[left + bankOffsetA];
    int saveRight = buffer[right + bankOffsetB];

    __syncthreads();
    for(int depth = size>>1; depth > 0; depth >>= 1){
        __syncthreads();
        if(threadId < depth){
            int forLeft = offset * (2 * threadId + 1) - 1;
            int forRight = offset * (2 * threadId + 2) - 1;
            forLeft += CONFLICT_FREE_OFFSET(forLeft);
            forRight += CONFLICT_FREE_OFFSET(forRight); 
            
            buffer[forRight] += buffer[forLeft];
        }
        offset *= 2;
    }

    if(threadId == 0) buffer[(size - 1) + CONFLICT_FREE_OFFSET(size - 1)] = 0;

    for(int depth = 1; depth < size; depth *= 2){
        offset /= 2;
        __syncthreads();
        if(threadId < depth){
            int forLeft = offset * (2 * threadId + 1) - 1;
            int forRight = offset*(2 * threadId + 2) - 1;
            forLeft += CONFLICT_FREE_OFFSET(forLeft);
            forRight += CONFLICT_FREE_OFFSET(forRight); 
            
            int tmp = buffer[forLeft];
            buffer[forLeft] = buffer[forRight];
            buffer[forRight] += tmp;
        }
    }

    __syncthreads();
    if(threadId < size / 2){
        histogram[left + blockOffset] = buffer[left + bankOffsetA] + saveLeft;
        histogram[right + blockOffset] = buffer[right + bankOffsetB] + saveRight;
        if(blockSum != NULL && threadIdx.x == (size / 2) - 1)
            blockSum[blockIdx.x] = buffer[right + bankOffsetB] + saveRight;
    }
}

__global__ void spreadSum(int *histogram, int *blockSum, int size){
    extern __shared__ int buffer[];
    int blockOffset = blockOffset = (blockIdx.x + 1) * size;

    if(threadIdx.x == 0)
        memcpy(buffer, histogram + blockOffset, size * sizeof(int));
    __syncthreads();

    int blockValue = blockSum[blockIdx.x];

    for(int i = threadIdx.x; i < size; i += blockDim.x)
        buffer[i] += blockValue;
    __syncthreads();

    if(threadIdx.x == 0)
        memcpy(histogram + blockOffset, buffer, size * sizeof(int));
    __syncthreads();
}

__global__ void countSort(int *array, int *histogram, int size){
    int taskPerBlock = (size > gridDim.x) ? (size / gridDim.x) : 1;
    int taskPerThread = (taskPerBlock > blockDim.x) ? (taskPerBlock / blockDim.x) : 1;
    int threadOffset = threadIdx.x + blockIdx.x * blockDim.x;
    int threadStart = threadOffset * taskPerThread;
    
    for(int i = threadStart; i < threadStart + taskPerThread; ++i){
        int j = (i == 0) ? 0 : histogram[i - 1];
        int k = histogram[i];
        for(; j < k; ++j){
            array[j] = i;
        }
    }

}

__host__ unsigned int maxTwo(int n){
    unsigned int x = n;
    x = x - 1;
    for(int i = 1; i < 32; i <<= 1){
        x |= x >> i;
    }
    return x + 1;
}