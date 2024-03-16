#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#define N 32 // blockDim.x = blockDim.y = N 
#define X (threadIdx.x) // to simplify expressions 
#define Y (threadIdx.y) // to simplify expressions 
#define RANGE_OF_SEARCH_SHRINK 0.5f
#define MAX_N_GENES 10000

#pragma pack(push, 1)
struct traces {
    int n_traces;
    int trace_lengths[MAX_N_GENES];
    float * trace_y[MAX_N_GENES]; // assume that the x values are integers from 0 to trace_lengths[i] - 1
    // for trace i, the coordiantes of j-th point is (j, trace_y[i][j])
};
#pragma pack(pop)

bool import(const char * const path, traces * imported, std::vector<std::string> & gene_names) {
    gene_names.clear();

    std::ifstream file_in(path);
    if (!file_in) {
        printf("\tFile '%s' doesn't exist.\n", path);
        return false;
    } 

    // Read lines from the file
    std::string line;
    int traceIdx = 0;
    while (std::getline(file_in, line)) {
        // Split the line by comma and convert to float
        std::istringstream iss(line);
        std::string token;

        // Skip first and last elements (assuming they are not needed)
        std::getline(iss, token, ','); // Skip first element
        gene_names.emplace_back(token);

        std::vector<float> values;
        while (std::getline(iss, token, ',')) {
            values.push_back(std::stof(token));
        }
        (imported->trace_lengths)[traceIdx] = values.size();

        (imported->trace_y)[traceIdx] = new float[imported->trace_lengths[traceIdx]];
        for(int i = 0; i < values.size(); ++i) {
            (imported->trace_y)[traceIdx][i] = values[i];
        }

        ++traceIdx;
    }
    
    imported->n_traces = traceIdx;

    std::cout << "\t" << gene_names.size() << " gene(s) imported." << std::endl;
    return true;
}

// placeholder
__device__ float score(float x, float y, float z) {
    // use cross-correlation as sceo
    return (float)(100.0f - (x - N/2 - z) * (x - N/2 - z) - (y - N/2 + z) * (y - N/2 + z));
}

__global__ void kernel(const traces *lib_traces, const float * exp, int exp_len, float * best_corr, float * best_offset, float * best_scaling) {
    // resolution of both theory and experiment is 1 nm per point. So, the total extension = data_size
    __shared__ float cache_correlation[N][N];
    __shared__ int   cache_offset_idx[N][N];
    __shared__ int   cache_scaling_idx[N][N];

    int lib_len = (lib_traces->trace_lengths)[blockIdx.x];
    float* lib_trace = (lib_traces->trace_y)[blockIdx.x];

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
        cache_correlation[X][Y] = score((float)X, (float)Y, blockIdx.x); 
        __syncthreads(); // wait for all correlation is calculated
        
        // My note for O(log n) redution of 2d array: because this is a 2D array, I have to compare 4 values each time in 4 smaller submatrices. 
        // Find max of each row has a time complexity of O(nlogn), therefore not right. 
        // I also need 2 arrays to track the argx and argy of max correlation. 
        float max1, xmax1, ymax1, max2, xmax2, ymax2;
        int i = blockDim.x/2;
        while (i != 0) {
            if (X < i && Y < i) {
                // printf("Kernel %d: X=%d,Y=%d,i=%d|x1=%f,x2=%f,x3=%f,x4=%f|X=%d,Y=%d,R=%f.\n",blockIdx.x, X, Y, i,
                // cache_correlation[X][Y], cache_correlation[X+i][Y],
                // cache_correlation[X][Y+i],cache_correlation[X + i][Y + i], 
                // cache_offset_idx [X][Y], cache_scaling_idx[X][Y], cache_correlation[X][Y]);
                max1 =  cache_correlation[X][Y    ] > cache_correlation[X + i][Y    ] ? cache_correlation[X][Y    ] : cache_correlation[X + i][Y    ];
                xmax1 = cache_correlation[X][Y    ] > cache_correlation[X + i][Y    ] ? cache_offset_idx [X][Y    ] : cache_offset_idx [X + i][Y    ];
                ymax1 = cache_correlation[X][Y    ] > cache_correlation[X + i][Y    ] ? cache_scaling_idx[X][Y    ] : cache_scaling_idx[X + i][Y    ];
                max2 =  cache_correlation[X][Y + i] > cache_correlation[X + i][Y + i] ? cache_correlation[X][Y + i] : cache_correlation[X + i][Y + i];
                xmax2 = cache_correlation[X][Y + i] > cache_correlation[X + i][Y + i] ? cache_offset_idx [X][Y + i] : cache_offset_idx [X + i][Y + i];
                ymax2 = cache_correlation[X][Y + i] > cache_correlation[X + i][Y + i] ? cache_scaling_idx[X][Y + i] : cache_scaling_idx[X + i][Y + i];

                cache_correlation[X][Y] = max1 > max2 ? max1 : max2;
                cache_offset_idx [X][Y] = max1 > max2 ? xmax1: xmax2;
                cache_scaling_idx[X][Y] = max1 > max2 ? ymax1: ymax2;
                //printf("Kernel %d: X=%d,Y=%d,R=%f.\n", blockIdx.x, cache_offset_idx[X][Y], cache_scaling_idx[X][Y], cache_correlation[X][Y]);
            }
            __syncthreads();
            i >>= 1;
        }

        // adjust the offset and scaling search range by half and repeat find the global max correlation
        max_offset  = cache_offset_idx[0][0]   + RANGE_OF_SEARCH_SHRINK * 0.5 * (max_offset  - min_offset);
        min_offset  = cache_offset_idx[0][0]   - RANGE_OF_SEARCH_SHRINK * 0.5 * (max_offset  - min_offset);
        max_scaling = cache_scaling_idx [0][0] + RANGE_OF_SEARCH_SHRINK * 0.5 * (max_scaling -  min_scaling);
        min_scaling = cache_scaling_idx [0][0] - RANGE_OF_SEARCH_SHRINK * 0.5 * (max_scaling -  min_scaling);
        
        break;
    }


    best_corr   [blockIdx.x] = cache_correlation[0][0];
    best_offset [blockIdx.x] = (float) cache_offset_idx[0][0]/(N - 1) * (max_offset - min_offset) + min_offset;
    best_scaling[blockIdx.x] = (float) cache_scaling_idx [0][0]/(N - 1) * (max_scaling -  min_scaling)  + min_scaling;


        
    if(X == 0 && Y == 0) {
        printf("\t'kernel' print: block %d, corr[%d][%d] = %f.\n", blockIdx.x, cache_offset_idx[0][0], cache_scaling_idx[0][0], cache_correlation[0][0]);
    }
}

#ifdef DLL
// the output function
extern "C" __declspec(dllexport) int cuda_main(const traces* lib_traces, const float * experiment, int exp_len) {
#else
int cuda_main(const traces* lib_traces, const float * experiment, int exp_len) {
#endif
    
    dim3 blockDim(N, N);
    int grid_size = lib_traces->n_traces;

    float * d_exp;
    cudaMalloc((void**)&d_exp, sizeof(float) * exp_len);
    cudaMemcpy(d_exp, experiment, sizeof(float) * exp_len, cudaMemcpyHostToDevice);

    traces* d_lib;
    cudaMalloc((void**)&d_lib, sizeof(traces));
    // cudaMemcpy(d_lib, lib_traces, sizeof(traces), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_lib->n_traces), &(lib_traces->n_traces), sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(&(d_lib->trace_lengths), &(lib_traces->trace_lengths), sizeof(int) * lib_traces->n_traces, cudaMemcpyHostToDevice);
    for(int i = 0; i < lib_traces->n_traces; ++i) {
        float * d_ptr;
        float * ptr = (lib_traces->trace_y)[i];
        int len = (lib_traces->trace_lengths)[i];
        cudaMalloc((void**)&d_ptr, sizeof(float) * len);
        printf("Hello World! len=%d.\n",len);
        cudaMemcpy(d_ptr, ptr, sizeof(float) * len, cudaMemcpyHostToDevice);????why this part doesn't work!!
        d_lib->trace_y[i] = d_ptr;
    }
    

    float * d_best_corr;
    float * d_best_offset;
    float * d_best_scaling;
    cudaMalloc((void**)&d_best_corr, sizeof(float) * grid_size);
    cudaMalloc((void**)&d_best_offset, sizeof(float) * grid_size);
    cudaMalloc((void**)&d_best_scaling, sizeof(float) * grid_size);

    printf("grid: %d.\n", grid_size);
    kernel<<<grid_size, blockDim>>>(d_lib, d_exp, exp_len, d_best_corr, d_best_offset, d_best_scaling);


    
    float * best_corr   = new float[grid_size];
    float * best_offset = new float[grid_size];
    float * best_scaling  = new float[grid_size];
    cudaMemcpy(best_corr,   d_best_corr,   sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_offset, d_best_offset, sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_scaling,  d_best_scaling,  sizeof(float) * grid_size, cudaMemcpyDeviceToHost);

    cudaFree(d_best_corr);
    cudaFree(d_best_offset);
    cudaFree(d_best_scaling);
    for(int i = 0; i < lib_traces->n_traces; ++i) {
        cudaFree(d_lib->trace_y[i]);
    }
    cudaFree(d_lib);
    cudaFree(d_exp);

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

#define TEST
#ifdef TEST
int main(int argc, char * argv[]) {
    //about();
    
    traces imported;
    std::vector<std::string> gene_names;
    import("example/H5alpha.txt", &imported, gene_names);
    cuda_main(&imported, (imported.trace_y)[2], (imported.trace_lengths)[2]);

    // destruct
    for(int i = 0; i < imported.n_traces; ++i) {
        delete[] imported.trace_y[i];
    }
    return 0;
}
#endif