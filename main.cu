#include <iostream>

__global__ void kernel(int int * corr, ) {
    // map from threadIdx/BlockIdx to pixel position
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    int z = threadIdx.y + blockIdx.y * blockDim.y;
    if (j >= max_j || z >= max_z) {
        return;
    }
    int offset = j + z * blockDim.x * gridDim.x;

int main() {

    return 0;
}