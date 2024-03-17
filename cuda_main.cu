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

// struct to store a set of the unzipping traces, can be experimental data or the lib of theoretical unzipping traces.
#pragma pack(push, 1)
struct traces {
    // Copying of this struct is very tricky!!!
    int n_traces;
    int *trace_lengths;
    float ** trace_y; // assume that the x values are integers from 0 to trace_lengths[i] - 1
    // for trace i, the coordiantes of j-th point is (j, trace_y[i][j])
};
#pragma pack(pop)
/*
// A general formula to write code to copy a float *** array:
// Theoreticaly the number of nested loops can continue when the dimension become larger and larger. The #loop is one less than the #asterisk
//ptr1, ptr2, ptr3 are all helper variables, they need to be released!

//float data[i][j][k]; // assume that this is the data and data dimensions

float *** ptr1 = new[i];
for(i){
    float ** ptr2 = new[j];
    for(j){
        float * ptr3;
        cudaMalloc(ptr3);
        cudaMemcpy(ptr3, data[i][j], length_of_dim_k); 
        ptr2[j] = ptr3;
    }
    cudaMalloc(ptr2);
    cudaMemcpy(ptr2, data[i], length_of_dim_j); 
    ptr1[i] = ptr2;
    delete[] ptr2;
}
float *** d_data;
cudaMalloc (d_data);
cudaMemcpy(d_data, ptr1, length_of_dim_i)
delete[] ptr1;

//d_data is the pointer to the array on GPU

*/
traces* copy_to_dev(const traces* h_t) {
     // very tedious somehow!
    int* helper1;
    cudaMalloc((void**)&helper1, h_t->n_traces * sizeof(int));
    cudaMemcpy(helper1, h_t->trace_lengths, h_t->n_traces * sizeof(int), cudaMemcpyHostToDevice);

    float ** helper2 = new float*[h_t->n_traces];// helper 2 is a ptr points to an array on host, the array is addresses to arrays which is also on device 
    for (int i = 0; i < h_t->n_traces; ++i) {
        // the 2 lines of code are correct, they are from my optimization of testing code. However, it is too brain xxxxx so I replaced it with clearer version!
        // cudaMalloc((void**)(helper2 + i), h_t->trace_lengths[i] * sizeof(float));
        // cudaMemcpy(helper2[i], (h_t->trace_y)[i], h_t->trace_lengths[i] * sizeof(float), cudaMemcpyHostToDevice);
        float * ptr;
        cudaMalloc((void**)&ptr, h_t->trace_lengths[i] * sizeof(float));
        cudaMemcpy(ptr, (h_t->trace_y)[i], h_t->trace_lengths[i] * sizeof(float), cudaMemcpyHostToDevice);
        (helper2[i]) = ptr;
        if (i == 1) {
            std::cout << "!!!!!!!!!!!!" << h_t->trace_y[i][2] << std::endl;
        }
    }

    float ** helper3;// helper 2 is a ptr points to an array on device, the array is addresses to arrays which is also on device
    cudaMalloc((void**)&helper3, h_t->n_traces * sizeof(float*));
    cudaMemcpy(helper3, helper2, h_t->n_traces * sizeof(float*), cudaMemcpyHostToDevice);
    delete[] helper2;

    // copy the value of "h_t->n_traces, helper1, helper3" to device, and save address in d_t
    traces h_temp {h_t->n_traces, helper1, helper3};
    traces* d_t;
    cudaMalloc((void**)&d_t, sizeof(traces));
    cudaMemcpy(d_t, &h_temp, sizeof(traces), cudaMemcpyHostToDevice);
    printf("Hello World Final!\n");
    return d_t;
}
void dev_del(traces* d_t) {
    for (int i = 0; i < d_t->n_traces; ++i){
        cudaFree(d_t->trace_y[i]);
    }
    cudaFree(d_t->trace_lengths);
    cudaFree(d_t->trace_y);
}
void host_del(traces * h_t){
    for (int i = 0; i < h_t->n_traces; ++i){
        delete[] h_t->trace_y[i];
    }
    delete[] h_t->trace_lengths;
    delete[] h_t->trace_y;
}

// wrapper to read unzipping data then create a traces object
bool import(const char * const path, traces * imported, std::vector<std::string> & gene_names, int start = 0, int end = 10000000000) {
    gene_names.clear();

    std::ifstream file_in(path);
    if (!file_in) {
        printf("\tFile '%s' doesn't exist.\n", path);
        return false;
    } 

    // Read lines from the file
    std::string line;

    int temp = 0;
    
    int *trace_lengths_temp = new int [100000];
    float ** trace_y_temp   = new float *[100000];

    while(++temp  <= start && std::getline(file_in, line)) {}
    int traceIdx = 0;
    while (traceIdx < (end - start) && std::getline(file_in, line)) {
        // Split the line by comma and convert to float
        std::istringstream iss(line);

        std::string token;

        // Skip first and last elements (assuming they are not needed)
        std::getline(iss, token, ','); // Skip first element
        gene_names.emplace_back(token);
        //std::cout << token << std::endl;

        std::vector<float> values;
        while (std::getline(iss, token, ',')) {
            values.push_back(std::stof(token));
        }
        (trace_lengths_temp)[traceIdx] = values.size();

        (trace_y_temp)[traceIdx] = new float[trace_lengths_temp[traceIdx]];
        for(int i = 0; i < values.size(); ++i) {
            trace_y_temp[traceIdx][i] = values[i];
        }

        ++traceIdx;
    }
    
    imported->n_traces = traceIdx;
    imported->trace_lengths = new int [traceIdx];
    imported->trace_y = new float * [traceIdx];
    memcpy(imported->trace_lengths, trace_lengths_temp, sizeof(int) * traceIdx);
    memcpy(imported->trace_y, trace_y_temp, sizeof(float*) * traceIdx);

    std::cout << "\t" << gene_names.size() << " gene(s) imported." << std::endl;
    for(int i = 0; i < imported->n_traces; ++i){
        std::cout << "\tgene: " << gene_names[i];
        printf("length = %d, trace_y[%d][2] = %f.\n", imported->trace_lengths[i], i, imported->trace_y[i][2]);
    }
    delete[] trace_lengths_temp;
    delete[] trace_y_temp;
    return true;
}

// placeholder
__device__ float score2(float x, float y, float z) {
    // use cross-correlation as sceo
    return (float)(100.0f - (x - N/2 - z) * (x - N/2 - z) - (y - N/2 + z) * (y - N/2 + z));
}
__device__ float score(float* exp, float* lib, int len_exp, int len_lib) {
    // use cross-correlation as score
    return 0;
}

__global__ void kernel(const traces *lib_traces, const traces *exp_traces, float * best_corr, float * best_offset, float * best_scaling) {
    // resolution of both theory and experiment is 1 nm per point. So, the total extension = data_size
    __shared__ float cache_correlation[N][N];
    __shared__ int   cache_offset_idx[N][N];
    __shared__ int   cache_scaling_idx[N][N];

    if (X == 0 && Y == 0 && blockIdx.x == 0 && blockIdx.y == 0) {
        printf("lib_traces->n_traces: %d.\n", lib_traces->n_traces);
        printf("lib_traces->trace_lengths[1]: %d.\n", lib_traces->trace_lengths[1]);
        printf("lib_traces->trace_y[1][2]: %f.\n", lib_traces->trace_y[1][2]);
    }

    int lib_len = (lib_traces->trace_lengths)[blockIdx.y];
    float* lib_trace = (lib_traces->trace_y)[blockIdx.y];
    int exp_len = (exp_traces->trace_lengths)[blockIdx.x];
    float* exp_trace = (exp_traces->trace_y)[blockIdx.x];

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
        cache_correlation[X][Y] = score2((float)X, (float)Y, blockIdx.x); 
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


    best_corr   [gridDim.x * blockIdx.y + blockIdx.x] = cache_correlation[0][0];
    best_offset [gridDim.x * blockIdx.y + blockIdx.x] = (float) cache_offset_idx[0][0]/(N - 1) * (max_offset - min_offset) + min_offset;
    best_scaling[gridDim.x * blockIdx.y + blockIdx.x] = (float) cache_scaling_idx [0][0]/(N - 1) * (max_scaling -  min_scaling)  + min_scaling;


    best_corr   [gridDim.x * blockIdx.y + blockIdx.x] = lib_trace[0];
    best_offset [gridDim.x * blockIdx.y + blockIdx.x] = exp_trace[0];
    best_scaling[gridDim.x * blockIdx.y + blockIdx.x] = lib_trace[0] + exp_trace[0];
        
    if(X == 0 && Y == 0) {
        //printf("\t'kernel' print: block %d,%d, corr[%d][%d] = %f.\n", blockIdx.x, blockIdx.y, cache_offset_idx[0][0], cache_scaling_idx[0][0], cache_correlation[0][0]);
        printf("\t'kernel' block (%d,%d)(%d), corr: %f, offset: %f, scaling: %f.\n", blockIdx.x, blockIdx.y, gridDim.x * blockIdx.y + blockIdx.x,
                                    best_corr   [gridDim.x * blockIdx.y + blockIdx.x], 
                                    best_offset [gridDim.x * blockIdx.y + blockIdx.x], 
                                    best_scaling[gridDim.x * blockIdx.y + blockIdx.x]);
    }
}

#ifdef DLL
// the output function
extern "C" __declspec(dllexport) int cuda_main(const traces* lib_traces, const float * experiment, int exp_len) {
#else
int cuda_main(const traces* lib_traces, const traces* exp_traces) {
#endif
    
    dim3 blockDim(N, N);
    dim3 gridDim(exp_traces->n_traces, lib_traces->n_traces);
    int grid_size = exp_traces->n_traces * lib_traces->n_traces;

    
    traces * d_exp = copy_to_dev(exp_traces);
    traces * d_lib = copy_to_dev(lib_traces);
    

    float * d_best_corr;
    float * d_best_offset;
    float * d_best_scaling;
    cudaMalloc((void**)&d_best_corr, sizeof(float) * grid_size);
    cudaMalloc((void**)&d_best_offset, sizeof(float) * grid_size);
    cudaMalloc((void**)&d_best_scaling, sizeof(float) * grid_size);

    printf("grid: %d.\n", grid_size);
    kernel<<<gridDim, blockDim>>>(d_lib, d_exp, d_best_corr, d_best_offset, d_best_scaling);


    
    float * best_corr   = new float[grid_size];
    float * best_offset = new float[grid_size];
    float * best_scaling  = new float[grid_size];
    cudaMemcpy(best_corr,   d_best_corr,   sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_offset, d_best_offset, sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_scaling,  d_best_scaling,  sizeof(float) * grid_size, cudaMemcpyDeviceToHost);

    // for (int i = 0; i < grid_size; ++i) {
    //     std::cout << ">>" << best_corr[i] << ',' << best_offset[i] << std::endl;
    // }

    cudaFree(d_best_corr);
    cudaFree(d_best_offset);
    cudaFree(d_best_scaling);
    dev_del(d_lib);
    dev_del(d_exp);

    return 0;
}

int about(){
    int device;
    cudaGetDevice(&device);

    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, device);

    printf("Maximum threads per block: %d\n", props.maxThreadsPerBlock);
    printf("Maximum dimensions of threads: %d x %d x %d\n", 
           props.maxThreadsDim[0], props.maxThreadsDim[1], props.maxThreadsDim[2]);
    printf("Maximum dimensions of grid: %d x %d x %d\n", 
           props.maxGridSize[0], props.maxGridSize[1], props.maxGridSize[2]);

    return 0;
}

#define TEST
#ifdef TEST
int main(int argc, char * argv[]) {
    about();
    
    traces mimic_exp {};
    std::vector<std::string> xxx;
    import("example/H5alpha.txt", &mimic_exp, xxx, 3, 6);

    traces library {};
    // std::vector<std::string> gene_names;
    // import("example/H5alpha.txt", &library, gene_names);


    cuda_main(&mimic_exp, &mimic_exp);

    // destruct
    host_del(&library);
    host_del(&mimic_exp);
    return 0;
}
#endif