#include <iostream>
#include <cuda_runtime.h>

// Define your struct
struct MyStruct {
    int size;       // Size of the array
    float* data;    // Pointer to float array
};

// Function to allocate memory for struct and data on the device
void allocateMemoryOnDevice(MyStruct* d_struct, MyStruct* h_struct) {
    // Copy struct to device
    std::cout << "check1!" << std::endl;
    // cudaMemcpy(d_struct, h_struct, sizeof(MyStruct), cudaMemcpyHostToDevice);

    std::cout << "check2!" << std::endl;
    // Allocate memory for data on device
    // cudaMalloc(&(d_struct->data), h_struct->size * sizeof(float));

    float* d_ptr;
    
    cudaMalloc((void**)&d_ptr, h_struct->size * sizeof(float));
    cudaMemcpy(d_ptr, h_struct->data, h_struct->size * sizeof(float), cudaMemcpyHostToDevice);
    std::cout << "check3!" << std::endl;
    // Copy data to device
    MyStruct *temp = new MyStruct;
    temp->size = h_struct->size;
    temp->data = d_ptr;
    cudaMemcpy(d_struct, temp, sizeof(MyStruct), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_struct->data, h_struct->data, h_struct->size * sizeof(float), cudaMemcpyHostToDevice);// I cannot dereference data! this won't work!
    std::cout << "check4!" << std::endl;
}

int main() {
    MyStruct h_struct;
    MyStruct *d_struct;

    std::cout << "start!" << std::endl;
    // Initialize h_struct (example)
    h_struct.size = 10;
    h_struct.data = new float[h_struct.size];
    for (int i = 0; i < h_struct.size; ++i) {
        h_struct.data[i] = static_cast<float>(i);
    }

    std::cout << "check!" << std::endl;
    // Allocate memory for struct on the device
    cudaMalloc(&d_struct, sizeof(MyStruct));

    std::cout << "check!" << std::endl;
    // Allocate memory for struct and data on the device
    allocateMemoryOnDevice(d_struct, &h_struct);
    std::cout << "check!" << std::endl;

    // Clean up
    delete[] h_struct.data;
    cudaFree(d_struct->data);
    cudaFree(d_struct);
    return 0;
}