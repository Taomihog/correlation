#include <iostream>
#include <cuda_runtime.h>
#include "cuda_main.h"

#define N 32 // blockDim.x = blockDim.y = N 
#define X (threadIdx.x) // to simplify expressions 
#define Y (threadIdx.y) // to simplify expressions 
#define RANGE_OF_SEARCH_SHRINK 0.5f

__global__ void correlation(){
    
}

__global__ void kernel(float * best_corr, float * best_offset, float * best_scaling, const float * const experiment, const float * const theory_lib, int data_size) {
    // resolution of both theory and experiment is 1 nm per point. So, the total extension = data_size
    __shared__ float cache_correlation[N][N];
    __shared__ int   cache_offset_idx[N][N];
    __shared__ int   cache_scaling_idx[N][N];


    float min_scaling = 0.9, max_scaling = 1.1;
    float min_offset = -500, max_offset =500;
    
    // The surface has local minima so I cannot use general minimization method, for example downhill simplex.
    // I have to use a brute force search, after each loop I shrink the search area a little bit.
    while(true) {
        // reset cache for this thread, the values will be changed during reduction to same the max_correlations' threadIdx.x and threadIdx.y
        cache_offset_idx[X][Y] = X;
        cache_scaling_idx [X][Y] = Y;
        // set offset and scaling of this thread
        float offset  = (float)X/(N - 1) * (max_offset - min_offset) + min_offset;
        float scaling = (float)Y/(N - 1) * (max_scaling  - min_scaling)  + min_scaling;

        // Calculate the correlation:
        float act_x = experiment[0] + theory_lib[blockIdx.x * data_size + 0] + offset + scaling; //actual extension
        cache_correlation[X][Y] = (float)1000 - (X - N/2 + blockIdx.x) * (X - N/2 + blockIdx.x) - (Y - N/2, 2) * (Y - N/2, 2); // placeholder

        __syncthreads(); // wait for all correlation is calculated

        // My note for O(log n) redution of 2d array: because this is a 2D array, I have to compare 4 values each time in 4 smaller submatrices. 
        // Find max of each row has a time complexity of O(nlogn), therefore not right. 
        // I also need 2 arrays to track the argx and argy of max correlation. 
        float max1, xmax1, ymax1, max2, xmax2, ymax2;
        int i = blockDim.x/2;
        while (i != 0) {
            if (X < i && Y < i) {
                max1 =  cache_correlation[X][Y] > cache_correlation[X][Y + i] ? cache_correlation[X][Y] : cache_correlation[X][Y + i];
                xmax1 = cache_correlation[X][Y] > cache_correlation[X][Y + i] ? cache_offset_idx [X][Y] : cache_offset_idx [X][Y + i];
                ymax1 = cache_correlation[X][Y] > cache_correlation[X][Y + i] ? cache_scaling_idx  [X][Y] : cache_scaling_idx  [X][Y + i];

                max2 =  cache_correlation[X + i][Y] > cache_correlation[X + i][Y + i] ? cache_correlation[X + i][Y] : cache_correlation[X + i][Y + i];
                xmax2 = cache_correlation[X + i][Y] > cache_correlation[X + i][Y + i] ? cache_offset_idx [X + i][Y] : cache_offset_idx [X + i][Y + i];
                ymax2 = cache_correlation[X + i][Y] > cache_correlation[X + i][Y + i] ? cache_scaling_idx  [X + i][Y] : cache_scaling_idx  [X + i][Y + i];

                cache_correlation[X][Y] = max1 > max2 ? max1 : max2;
                cache_offset_idx [X][Y] = max1 > max2 ? xmax1: xmax2;
                cache_scaling_idx  [X][Y] = max1 > max2 ? ymax1: ymax2;
            }
            __syncthreads();
            i >>= 1;
        }
        
        break;

        // adjust the offset and scaling search range by half and repeat find the global max correlation
        max_offset = cache_offset_idx[0][0] + RANGE_OF_SEARCH_SHRINK * 0.5 * (max_offset - min_offset);
        min_offset = cache_offset_idx[0][0] - RANGE_OF_SEARCH_SHRINK * 0.5 * (max_offset - min_offset);
        max_scaling  = cache_scaling_idx [0][0] + RANGE_OF_SEARCH_SHRINK * 0.5 * (max_scaling -  min_scaling);
        min_scaling  = cache_scaling_idx [0][0] - RANGE_OF_SEARCH_SHRINK * 0.5 * (max_scaling -  min_scaling);
    }


    best_corr  [blockIdx.x] = cache_correlation[0][0];
    best_offset[blockIdx.x] = (float) cache_offset_idx[0][0]/(N - 1) * (max_offset - min_offset) + min_offset;
    best_scaling [blockIdx.x] = (float) cache_scaling_idx [0][0]/(N - 1) * (max_scaling -  min_scaling)  + min_scaling;

        
    if(X == 0 && Y == 0) {
        printf("In kernel print, block %d, corr[%d][%d] = %f.\n", blockIdx.x, cache_offset_idx[0][0], cache_scaling_idx[0][0], cache_correlation[0][0]);
    }
}

int cuda_main() {
    
    dim3 blockDim(16, 16);
    // placeholder;
    int Nblock = 10;//how many genes in the library
    int trace_length = 100;
    float * experiment = new float[trace_length]; 
    float * theory_lib = new float [Nblock * trace_length];
    // pad the empty data with -1.0f
    for (int i = 0; i < trace_length; ++i) {
        experiment[i] = -1.0f;
        for(int j = 0; j < Nblock; ++j) {
            theory_lib[i * Nblock + j] = -1.0f;
        }
    }

    float * d_experiment, * d_theory_lib;
    cudaMalloc((void**)&d_experiment, sizeof(float) * Nblock);
    cudaMalloc((void**)&d_theory_lib, sizeof(float) * Nblock * trace_length);
    cudaMemcpy(d_experiment, experiment, sizeof(float) * trace_length, cudaMemcpyHostToDevice);
    cudaMemcpy(d_theory_lib, theory_lib, sizeof(float) * Nblock * trace_length, cudaMemcpyHostToDevice);
    // end placeholder

    float * d_best_corr, * d_best_offset, * d_best_scaling;
    cudaMalloc((void**)&d_best_corr,   sizeof(float) * Nblock);
    cudaMalloc((void**)&d_best_offset, sizeof(float) * Nblock);
    cudaMalloc((void**)&d_best_scaling,  sizeof(float) * Nblock);

    

    // float * best_corr, float * best_scaling, float * best_offset, float * experiment, float ** theory_lib, int data_size
    kernel<<<Nblock, blockDim>>>(d_best_corr, d_best_offset, d_best_scaling, d_experiment, d_theory_lib, trace_length);


    
    float * best_corr   = new float[Nblock];
    float * best_offset = new float[Nblock];
    float * best_scaling  = new float[Nblock];
    cudaMemcpy(best_corr,   d_best_corr,   sizeof(float) * Nblock, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_offset, d_best_offset, sizeof(float) * Nblock, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_scaling,  d_best_scaling,  sizeof(float) * Nblock, cudaMemcpyDeviceToHost);

    cudaFree(d_best_corr);
    cudaFree(d_best_offset);
    cudaFree(d_best_scaling);
    cudaFree(d_experiment);
    cudaFree(d_theory_lib);

    delete[] experiment;
    delete[] theory_lib;
    return 0;
}

int about(){
    
    int device;
    cudaGetDevice(&device);

    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);

    printf("Maximum threads per block: %d\n", props.maxThreadsPerBlock);
    printf("Maximum dimensions of grid: %d x %d x %d\n", 
           props.maxGridSize[0], props.maxGridSize[1], props.maxGridSize[2]);

    return 0;
}

#ifdef TEST
int main() {
    about();
    cuda_main();
    return 0;
}
#endif