// Command to build: 
// nvcc main.cu -o main
// Run executable:
// main example/H5alpha.csv example/H5alpha.csv

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#define N 32 // blockDim.x = blockDim.y = N 
#define X (threadIdx.x) // to simplify expressions 
#define Y (threadIdx.y) // to simplify expressions 
#define RANGE_OF_SEARCH_SHRINK 0.3f
// #define TEST_PRINT
// #define TEST_PRINT2


// data struct and helper functions for copy and delete
#pragma pack(push, 1)
// struct to store a set of the unzipping traces, can be experimental data or the lib of theoretical unzipping traces.
struct traces {
    int n_traces;
    int *trace_lengths;
    float ** trace_y; // assume that the x values are integers from 0 to trace_lengths[i] - 1
    // Therefore, the coordiante of j-th point of the i-th trace's is (j, trace_y[i][j]), where i < n_traces and j < trace_lengths[n_traces].
};
#pragma pack(pop)
traces* copy_to_dev(const traces* h_t) {
/*
// Very tricky to copy a ptr array of ptrs (of ptrs of ptrs...)!
// Theoreticaly the #loop is one less than the #asterisk
//ptr1, ptr2, ptr3 are all helper variables, they need to be released!
// Here is an exmple of 3D array (float ***, #asterisk is 3, need a nested for loop of depth = 2 to copy). 
// I am not sure if there is better way to do it.

//float data[i][j][k]; // data
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
float *** d_data;//d_data is the pointer to the array on GPU
cudaMalloc (d_data);
cudaMemcpy(d_data, ptr1, length_of_dim_i)
delete[] ptr1;
*/
    int* helper1;
    cudaMalloc((void**)&helper1, h_t->n_traces * sizeof(int));
    cudaMemcpy(helper1, h_t->trace_lengths, h_t->n_traces * sizeof(int), cudaMemcpyHostToDevice);

    float ** helper2 = new float*[h_t->n_traces];
    for (int i = 0; i < h_t->n_traces; ++i) {float * ptr;
        cudaMalloc((void**)&ptr, h_t->trace_lengths[i] * sizeof(float));
        cudaMemcpy(ptr, (h_t->trace_y)[i], h_t->trace_lengths[i] * sizeof(float), cudaMemcpyHostToDevice);
        (helper2[i]) = ptr;
    }

    float ** helper3;
    cudaMalloc((void**)&helper3, h_t->n_traces * sizeof(float*));
    cudaMemcpy(helper3, helper2, h_t->n_traces * sizeof(float*), cudaMemcpyHostToDevice);
    delete[] helper2;

    traces h_temp {h_t->n_traces, helper1, helper3};
    traces* d_t;
    cudaMalloc((void**)&d_t, sizeof(traces));
    cudaMemcpy(d_t, &h_temp, sizeof(traces), cudaMemcpyHostToDevice);
    return d_t;
}
__device__ void dev_del(traces* d_t) {
    //don't use, it need to work in a "copy to dev" way.
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
bool import(const char * const path, traces * imported, std::vector<std::string> & gene_names, int start = 0, int end = 10000000) {
    gene_names.clear();

    std::ifstream file_in(path);
    if (!file_in) {
        printf("\tFile '%s' doesn't exist.\n", path);
        return false;
    } 

    int *trace_lengths_temp = new int [100000];
    float ** trace_y_temp   = new float *[100000];
    int temp = 0;
    std::string line;
    while(++temp  <= start && std::getline(file_in, line)) {}
    int traceIdx = 0;
    while (traceIdx < (end - start) && std::getline(file_in, line)) {
        std::istringstream iss(line);
        std::string token;

        std::getline(iss, token, ','); // The first element is trace name
        gene_names.emplace_back(token);

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

    #ifdef TEST_PRINT
    std::cout << "\tTEST_PRINT (import): " << std::endl;
    for(int i = 0; i < imported->n_traces; ++i){
        std::cout << "\tgene: " << gene_names[i];
        printf("length = %d, trace_y[%d][0] = %f.\n", imported->trace_lengths[i], i, imported->trace_y[i][0]);
    }
    #endif

    delete[] trace_lengths_temp;
    delete[] trace_y_temp;
    return true;
}

__device__ float score_func(const float* exp, const float* lib, int len_exp, int len_lib, float offset, float scaling) {
    // use a slightly modified cross-correlation as score function to align 2 traces
    float j_frac;
    int j = 0;
    int i = 0;

    float exp_i = 0.0f;
    float ac_lib = 0.0f;
    float ac_exp = 0.0f;
    float cc = 0.0f;

    for(; i < len_lib; ++i){
        j_frac = (i - offset) / scaling; // = (x_theory - offset) / scaling = x coordinate of experimental trace
        if (j_frac < 0) {
            continue;
        }
        if (j_frac >= len_exp) {
            break;
        }
        j = j_frac;
        exp_i = exp[j] * (j + 1 - j_frac) + exp[j + 1] * (j_frac - j); // interp experimental data to theory
        ac_lib += lib[i] * lib[i];
        ac_exp += exp_i * exp_i;
        cc += lib[i] * exp_i; //interplation
    }
    float score = cc * cc / ac_lib / ac_exp;
    float panalty = ((i + 1)/ len_lib) * ((j + 1)  / len_exp); // there should be a penalty when less points overlap.

    return score * (0.9f + 0.1f * panalty) * (2.0f / (scaling + 1.0f / scaling));
}

__global__ void zero_the_force(traces * ptr) {
    int len = (ptr->trace_lengths)[blockIdx.y];
    // float* lib_trace = (float*)malloc(sizeof(float) * lib_len);
    // memcpy(lib_trace,(lib_traces->trace_y)[blockIdx.y], sizeof(float) * lib_len);
    float* trace = (ptr->trace_y)[blockIdx.y];
    float avg = 0.0f;
    float y0 = trace[0];
    for(int i = 0; i < len; ++i){
        trace[i] -= y0;
        avg += trace[i];
    }
    avg /= len;
    for(int i = 0; i < len; ++i){
        trace[i] -= avg;
    }
}

//calculate the 
__global__ void kernel(const traces *lib_traces, const traces *exp_traces, float * best_corr, float * best_offset, float * best_scaling) {
    // resolution of both theory and experiment is 1 nm per point. So, the total extension = data_size
    __shared__ float cache_correlation[N][N];
    __shared__ int   cache_offset_idx[N][N];
    __shared__ int   cache_scaling_idx[N][N];
    __shared__ float min_scaling;
    __shared__ float max_scaling;
    __shared__ float min_offset;
    __shared__ float max_offset;


    int lib_len = (lib_traces->trace_lengths)[blockIdx.y];
    // float* lib_trace = (float*)malloc(sizeof(float) * lib_len);
    // memcpy(lib_trace,(lib_traces->trace_y)[blockIdx.y], sizeof(float) * lib_len);
    float* lib_trace = (lib_traces->trace_y)[blockIdx.y];

    int exp_len = (exp_traces->trace_lengths)[blockIdx.x];
    // float* exp_trace = (float*)malloc(sizeof(float) * exp_len);
    // memcpy(exp_trace,(exp_traces->trace_y)[blockIdx.y], sizeof(float) * exp_len);
    float* exp_trace = (exp_traces->trace_y)[blockIdx.x];
    
    if (X == 0 && Y == 0) {
        // not necessary to use __shared__
        // set scaling and offset
        min_scaling = 0.95f;
        max_scaling = 1.05f;
        min_offset  = -200.0f;
        max_offset  = 200.0f;
    }
    __syncthreads();
    
        
    #ifdef TEST_PRINT
    if(X == 0 && Y == 0){
        printf("\tblock[%d][%d]: lib_trace: %f, length: %d, ac_lib:%f. exp_trace: %f, length: %d. ac_exp:%f.\n", 
            blockIdx.x, blockIdx.y, lib_trace[0],lib_len, ac_lib, exp_trace[0],exp_len, ac_exp);
    }
    #endif

    // The surface has local minima so I cannot use general minimization method, for example downhill simplex.
    // I have to use a brute force search, after each loop I shrink the search area a little bit.
    const int max_iteration = 20;
    int iteration = 0;
    while(true) {
        // reset cache for this thread, the values will be changed during reduction to same the max_correlations' threadIdx.x and threadIdx.y
        cache_offset_idx[X][Y] = X;
        cache_scaling_idx [X][Y] = Y;
        // set offset and scaling of this thread
        float offset  = (float)X/(N - 1) * (max_offset - min_offset) + min_offset;
        float scaling = (float)Y/(N - 1) * (max_scaling  - min_scaling)  + min_scaling;

        // Calculate the correlation:
        cache_correlation[X][Y] = score_func(exp_trace, lib_trace, exp_len, lib_len, offset, scaling); 
        // cache_correlation[X][Y] = score_for_test((float)X, (float)Y, blockIdx.x); 
        __syncthreads(); // wait for all correlation is calculated
        
        // My note for O(log n) redution of 2d array: because this is a 2D array, I have to compare 4 values each time in 4 smaller submatrices. 
        // Find max of each row has a time complexity of O(nlogn), therefore not right. 
        // I also need 2 arrays to track the argx and argy of max correlation. 
        float max1, xmax1, ymax1, max2, xmax2, ymax2;
        int i = blockDim.x/2;
        while (i != 0) {
            if (X < i && Y < i) {
                // printf("Kernel %d: X=%d,Y=%d,i=%d|x1=%f,x2=%f,x3=%f,x4=%f|X=%d,Y=%d,R=%f.\n",blockIdx.x, X, Y, i,
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


        if (++iteration == max_iteration) {
            break; // must break before reducing offset and scaling, otherwise the calculation later will be incorrect
        }

        // adjust the offset and scaling search range by half and repeat find the global max correlation
        if(X == 0 && Y == 0) {
            float best_offset = (float)cache_offset_idx[0][0]/(N - 1) * (max_offset - min_offset) + min_offset;
            float new_offset_range = RANGE_OF_SEARCH_SHRINK * 0.5 * (max_offset  - min_offset);
            max_offset  = best_offset + new_offset_range;
            min_offset  = best_offset - new_offset_range;

            float best_scaling = (float)cache_scaling_idx [0][0]/(N - 1) * (max_scaling  - min_scaling)  + min_scaling;
            float new_scaling_range = RANGE_OF_SEARCH_SHRINK * 0.5 * (max_scaling -  min_scaling);
            max_scaling = best_scaling + new_scaling_range;
            min_scaling = best_scaling - new_scaling_range;

            // printf("max_offset: %f, min_offset: %f, max_scaling: %f, min_scaling: %f.\n", max_offset, min_offset, max_scaling, min_scaling);
        }
        __syncthreads();
    }

    best_corr   [gridDim.x * blockIdx.y + blockIdx.x] = cache_correlation[0][0];
    best_offset [gridDim.x * blockIdx.y + blockIdx.x] = ((float) cache_offset_idx[0][0]/(N - 1)) * (max_offset - min_offset) + min_offset;
    best_scaling[gridDim.x * blockIdx.y + blockIdx.x] = ((float) cache_scaling_idx [0][0]/(N - 1)) * (max_scaling -  min_scaling)  + min_scaling;
    

    #ifdef TEST_PRINT2
    if(X == 0 && Y == 0) {
        //printf("\t'kernel' print: block %d,%d, corr[%d][%d] = %f.\n", blockIdx.x, blockIdx.y, cache_offset_idx[0][0], cache_scaling_idx[0][0], cache_correlation[0][0]);
        printf("\t'kernel' block (%d,%d)(%d), corr: %f, offset: %f, scaling: %f.\n", blockIdx.x, blockIdx.y, gridDim.x * blockIdx.y + blockIdx.x,
                                    best_corr   [gridDim.x * blockIdx.y + blockIdx.x], 
                                    best_offset [gridDim.x * blockIdx.y + blockIdx.x], 
                                    best_scaling[gridDim.x * blockIdx.y + blockIdx.x]);
    }
    #endif

}

#ifdef DLL
extern "C" __declspec(dllexport) float*** align_traces(const traces* lib_traces, const float * experiment, int exp_len) {
#else
float*** align_traces(const traces* lib_traces, const traces* exp_traces) {
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



    zero_the_force<<<1, lib_traces->n_traces>>>(d_lib);
    zero_the_force<<<1, exp_traces->n_traces>>>(d_exp);
    kernel<<<gridDim, blockDim>>>(d_lib, d_exp, d_best_corr, d_best_offset, d_best_scaling);



    float * best_corr   = new float[grid_size];
    float * best_offset = new float[grid_size];
    float * best_scaling  = new float[grid_size];
    cudaMemcpy(best_corr,   d_best_corr,   sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_offset, d_best_offset, sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaMemcpy(best_scaling,  d_best_scaling,  sizeof(float) * grid_size, cudaMemcpyDeviceToHost);
    cudaFree(d_best_corr);
    cudaFree(d_best_offset);
    cudaFree(d_best_scaling);

    // don't work for now.
    // dev_del(d_lib);
    // dev_del(d_exp);
    
    float** best_corr_2d    = new float*[gridDim.x];
    float** best_offset_2d  = new float*[gridDim.x];
    float** best_scaling_2d = new float*[gridDim.x];
    for (int x = 0; x < gridDim.x; ++x) {
        best_corr_2d   [x] = new float[gridDim.y];
        best_offset_2d [x] = new float[gridDim.y];
        best_scaling_2d[x] = new float[gridDim.y];
        for (int y = 0; y < gridDim.y; ++y) {
            best_corr_2d   [x][y] = best_corr   [y * gridDim.x + x];
            best_offset_2d [x][y] = best_offset [y * gridDim.x + x];
            best_scaling_2d[x][y] = best_scaling[y * gridDim.x + x];
        }
    }
    delete[] best_corr;
    delete[] best_offset;
    delete[] best_scaling;
    
    float *** res = new float **[3];
    res[0] = best_corr_2d;
    res[1] = best_offset_2d;
    res[2] = best_scaling_2d;
    return res;
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

bool save_csv_then_free(const char* const filename_out, float** data, size_t rows, size_t cols) {
    
    std::ofstream file(filename_out);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file for writing." << std::endl;
        return false;
    }
    // Write data to CSV file
    for (size_t x = 0; x < rows; ++x) {
        float * row = data[x];
        for (size_t y = 0; y < cols; ++y) {
            file << row[y];
            if (y < cols - 1) {
                file << ",";
            }
        }
        delete [] row;
        file << std::endl;
    }
    delete[] data;
    file.close();
    return true;
}

#define TEST
#ifdef TEST
int main(int argc, const char * const argv[]) {
    // about();
    // const char * const temp = "example/H5alpha.txt";
    if (argc < 3) {
        printf("requires 2 file names: 1) experimental data and 2) theory trace lib");
        return 1;
    }

    traces exp_data {};
    std::vector<std::string> exp_trace_names;
    import(argv[1], &exp_data, exp_trace_names, 1, 10);

    traces theory_lib {};
    std::vector<std::string> theory_lib_gene_names;
    import(argv[2], &theory_lib, theory_lib_gene_names, 2, 11);
    for(auto gname:theory_lib_gene_names) {
        std::cout<< gname << std::endl;
    }



    float*** result = align_traces(&exp_data, &theory_lib);



    save_csv_then_free("corr.csv",result[0], exp_data.n_traces, theory_lib.n_traces);
    save_csv_then_free("offset.csv",result[1], exp_data.n_traces, theory_lib.n_traces);
    save_csv_then_free("scaling.csv",result[2], exp_data.n_traces, theory_lib.n_traces);

    delete[] result;
    host_del(&theory_lib);
    host_del(&exp_data);
    return 0;
}
#endif