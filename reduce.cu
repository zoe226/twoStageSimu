#include "cuda_runtime.h"
#include "cuda_runtime_api.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

#define BLOCK_SIZE 128

void random_init(int *data, int size);
void cpuTest(int& dout, int *din, int size);
__global__ void test_kernel(int *d_out, int *d_in);
__global__ void test_kernel1(int *d_out, int *d_in);
__global__ void test_kernel2(int *d_out, int *d_in);
__global__ void test_kernel3(int *d_out, int *d_in);
__global__ void test_kernel4(int *d_out, int *d_in);

int main()
{
    cudaError_t cudaError;
    cudaError = cudaSetDevice(7);
    if (cudaError == cudaSuccess)
    {
        std::cout<<"choose success"<<std::endl;   
    }
    else
    {
        std::cout<<"choose fail"<<std::endl;
    }

    const int arraySize = 16*1024*1024;
    int n_iter = 9;
    int *h_in;
    int h_out = 0;
    cudaMallocHost(&h_in,arraySize*sizeof(int));
    random_init(h_in,arraySize);

    double dur;
    clock_t start1,end1;
    start1 = clock();
    cpuTest(h_out,h_in,arraySize);
    end1 = clock();
    dur = (double)(end1 - start1);
    printf("cpu time: %f\n", (dur/CLOCKS_PER_SEC));

    int block_num = (arraySize + BLOCK_SIZE - 1) / BLOCK_SIZE / 2;
    dim3 blockDim(BLOCK_SIZE, 1, 1);
    dim3 gridDim(block_num, 1, 1);
    int block_num2 = block_num / BLOCK_SIZE / 2;
    dim3 gridDim2(block_num2, 1, 1);
    dim3 gridDim3(1, 1, 1);
    int *d_in, *d_out_L1, *d_out_L2, *d_out_L3;
    cudaMalloc(&d_in, arraySize * sizeof(int));
    cudaMalloc(&d_out_L1, block_num * sizeof(int));
    cudaMalloc(&d_out_L2, block_num2 * sizeof(int));
    cudaMalloc(&d_out_L3, 1 * sizeof(int));

    cudaMemcpy(d_in, h_in, arraySize * sizeof(int), cudaMemcpyDefault);

    test_kernel4 <<< gridDim, blockDim >>> (d_out_L1, d_in); // level one
    test_kernel4 <<< gridDim2, blockDim >>> (d_out_L2, d_out_L1); // level two
    test_kernel4 <<< gridDim3, blockDim >>> (d_out_L3, d_out_L2); // level three

    cudaEvent_t start,end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);

    cudaEventRecord(start);
    for (size_t i = 0; i < n_iter; i++)
    {
        test_kernel4 <<< gridDim, blockDim >>> (d_out_L1, d_in); // level one
        test_kernel4 <<< gridDim2, blockDim >>> (d_out_L2, d_out_L1); // level two
        test_kernel4 <<< gridDim3, blockDim >>> (d_out_L3, d_out_L2); // level three
    }
    cudaEventRecord(end);
    cudaEventSynchronize(end);
    flaot msec, sec;
    cudaEventElapsedTime(&msec, start, end);
    sec = msec / 1000.0 / n_iter;

    cudaEventDestroy(start);
    cudaEventDestroy(end);

    printf("Latency: %f\n", sec);
    
    int *h_d_out1, *h_d_out2, *h_d_out3;
    cudaMallocHost(&h_d_out1, block_num * sizeof(int));
    cudaMallocHost(&h_d_out2, block_num2 * sizeof(int));
    cudaMallocHost(&h_d_out3, sizeof(int));
    cudaMemcpy(h_d_out1, d_out_L1, block_num * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_d_out2, d_out_L2, block_num2 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_d_out3, d_out_L3, sizeof(int), cudaMemcpyDeviceToHost);

    int d_result1 = 0, d_result2 = 0;
    for (size_t i = 0; i < block_num; i++)
    {
        d_result1 += h_d_out1[i];
    }
    for (size_t i = 0; i < block_num2; i++)
    {
        d_result2 += h_d_out2[i];
    }

    cudaFreeHost(h_in);
    cudaFreeHost(h_d_out1);
    cudaFreeHost(h_d_out2);
    cudaFreeHost(h_d_out3);
    cudaFree(d_out_L1);
    cudaFree(d_out_L2);
    cudaFree(d_out_L3);
    
    return 0;
}

void random_init(int *data, int size)
{
    for (size_t i = 0; i < size; i++)
    {
        data[i] = int(rand() % (10 - 1)) + 1;
    }
}

void cpuTest(int& dout, int *din, int size)
{
    for(size_t i = 0; i < size; i++)
    {
        dout += din[i];
    }
}

__global__ void test_kernel(int *d_out, int *d_in)
{
    __shared__ int sdata[BLOCK_SIZE];

    unsigned int did = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    sdata[tid] = d_in[did];
    __syncthreads();

    for(unsigned int stride = 1; stride < blockDim.x; stride*=2)
    {
        if(tid % (2*stride) == 0)
        {
            sdata[tid] += sdata[tid + stride];
        }
        __syncthreads();
    }
    if(tid == 0)
    {
        d_out[blockIdx.x] = sdata[tid];
    }
}

__global__ void test_kernel1(int *d_out, int *d_in)
{
    __shared__ int sdata[BLOCK_SIZE];

    unsigned int did = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    sdata[tid] = d_in[did];
    __syncthreads();

    for(unsigned int stride = 1; stride < blockDim.x; stride*=2)
    {
        int index = stride * 2 * tid;
        if(index < blockDim.x)
        {
            sdata[index] += sdata[index + stride];
        }
        __syncthreads();
    }
    if(tid == 0)
    {
        d_out[blockIdx.x] = sdata[tid];
    }
}

__global__ void test_kernel2(int *d_out, int *d_in)
{
    __shared__ int sdata[BLOCK_SIZE];

    unsigned int did = blockDim.x * blockIdx.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    sdata[tid] = d_in[did];
    __syncthreads();

    for(unsigned int stride = blockDim.x/2; stride > 0; stride>>=1)
    {
        if(tid < stride)
        {
            sdata[tid] += sdata[tid + stride];
        }
        __syncthreads();
    }
    if(tid == 0)
    {
        d_out[blockIdx.x] = sdata[tid];
    }
}

__global__ void test_kernel3(int *d_out, int *d_in)
{
    __shared__ int sdata[BLOCK_SIZE];

    unsigned int did = (blockDim.x*2) * blockIdx.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    sdata[tid] = d_in[did] + d_in[did + blockDim.x];
    __syncthreads();

   for(unsigned int stride = blockDim.x/2; stride > 0; stride>>=1)
    {
        if(tid < stride)
        {
            sdata[tid] += sdata[tid + stride];
        }
        __syncthreads();
    }
    if(tid == 0)
    {
        d_out[blockIdx.x] = sdata[tid];
    }
}

__global__ void test_kernel4(int *d_out, int *d_in)
{
    __shared__ int sdata[BLOCK_SIZE];

    unsigned int did = (blockDim.x*2) * blockIdx.x + threadIdx.x;
    unsigned int tid = threadIdx.x;
    sdata[tid] = d_in[did] + d_in[did + blockDim.x];
    __syncthreads();

   for(unsigned int stride = blockDim.x/2; stride > 32; stride>>=1)
    {
        if(tid < stride)
        {
            sdata[tid] += sdata[tid + stride];
        }
        __syncthreads();
    }
    if(tid < 32)
    {
        sdata[tid] += sdata[tid + 32];
        sdata[tid] += sdata[tid + 16];
        sdata[tid] += sdata[tid + 8];
        sdata[tid] += sdata[tid + 4];
        sdata[tid] += sdata[tid + 2];
        sdata[tid] += sdata[tid + 1];
    }
    if(tid == 0)
    {
        d_out[blockIdx.x] = sdata[tid];
    }
}