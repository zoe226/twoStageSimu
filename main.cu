#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "device_functions.h"
#include "cuda_runtime_api.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <iostream>
#include <time.h>

void random_init(float *data, size_t size);
void cpuSgemm(float *a, float *b, float *c, const int M, const int N, const int K);
__global__ void naiveSgemm(float *a, float *b, float *c, const int M, const int N, const int K);
__global__ void sgemm_kernel(float *A, float *B, float *C, const int M, const int N, const int K);

int main()
{
	int m = 512;
	int n = 512;
	int k = 512;
	int n_iter = 9;

	cudaError_t cudaError;
	cudaError = cudaSetDevice(7);
	if(cudaError == cudaSuccess)
	{
		std::cout << "choose success" << std::endl;
	}
	else
	{
		std::cout << "choose fail" << std::endl;
	}
	int deviceId; 
	cudaError = cudaGetDevice(&deviceId);

	float *h_A, *h_B, *h_C, *h_d_C;
	cudaMallocHost(&h_A, m*k*sizeof(float));
	cudaMallocHost(&h_B, k*n*sizeof(float));
	cudaMallocHost(&h_C, m*n*sizeof(float));
	cudaMallocHost(&h_d_C, m*n*sizeof(float));
	random_init(h_A, m*k);
	random_init(h_B, k*n);

	double dur;
	clock_t start1,endl;
	start1 = clock();
	cpuSgemm(h_A, h_B, h_C, m, n, k);
	end1 = clock();
	dur = (double)(end1 - start1);
	printf("cpu time: %f\n",(dur/CLOCLS_PER_SEC));

	float *d_A, *d_B, *d_C;
	cudaMalloc(&d_A, m*k*sizeof(float));
	cudaMalloc(&d_B, k*n*sizeof(float));
	cudaMalloc(&d_C, m*n*sizeof(float));

	cudaMemcpy(d_A, h_A, m*k*sizeof(float),cudaMemcpyDefault);
	cudaMemcpy(d_B, h_B, k*n*sizeof(float),cudaMemcpyDefault);
	
	dim3 blockDim(32,32);
	dim3 gridDime((n+32-1)/32,(m+32-1)/32);

	naiveSgemm<<<gridDim,blockDim>>>(d_A,d_B,d_C,m,n,k);
	sgemm_kernel<<<gridDim,blockDim>>>(d_A,d_B,d_C,m,n,k);

	cudaEvent_t start,end;
	cudaEventCreate(&start);
	cudaEventCreate(&end);
	cudaEventRecord(start);

	for(size_t i = 0; i < n_iter; i++)
	{
		naiveSgemm<<<gridDim,blockDim>>>(d_A,d_B,d_C,m,n,k);
		sgemm_kernel<<<gridDim,blockDim>>>(d_A,d_B,d_C,m,n,k);
	}

	cudaEventRecord(end);
	cudaEventSynchronize(end);
	float msec, sec;
	cudaEventElapsedTime(&msec, start, end);
	sec = msec / 1000.0 / n_iter;

	cudaEventDestroy(start);
	cudaEventDestroy(end);

	printf("Latency: %f\n",sec);

	cudaMemcpy(h_d_C, d_C, m*n*sizeof(float),cudaMemcpyDeviceToHost);
	float abserror = 0.0;
	for(size_t i = 0; i < m*n; i++)
	{
		float temperror = abs(h_d_C[i]-h_C[i]);
		if(temperror > abserror)
		{
			abserror  = temperror;
		}
	}

	cudaFree(d_A);
	cudaFree(d_B);
	cudaFree(d_C);
	cudaFreeHost(h_A);
	cudaFreeHost(h_B);
	cudaFreeHost(h_C);
	cudaFreeHost(h_d_C);

	return 0;
}

void randomA_init(float *data, size_t size)
{
	for(size_t i = 0; i < size; i++)
	{
		data[i] = float(rand()) / RAND_MAX;
	}
}

void cpuSgemm(float *a, float *b, float *c, const int M, const int N, const int K)
{
	for(size_t m = 0; m < M; m++)
	{
		for(size_t n = 0; n < N; n++)
		{
			float psum = 0.0;
			for(size_t k = 0; k < K; k++)
			{
				psum += a[k + m*K] * b[n + k*N];
			}
			c[n+m*N] = psum;
		}
	}
}

__global__ void naiveSgemm(float *a, float *b, float *c, const int M, const int N, const int K)
{
	int n = blockIdx.x * blockDim.x + threadIdx.x;
	int m = blockIdx.y * blockDim.y + threadIdx.y;
	if(m < M && n < N)
	{
		float psum = 0.0;
		for(size_t k = 0; k < K; k++)
		{
			psum += a[k + m*K] * b[n+k*N];
		}
		c[n + m*N] = psum;
	}
}

__global__ void sgemm_kernel(float *A, float *B, float *C, const int M, const int N, const int K)
{
	// blocked matrix multiply
	__shared__ float tileA[32][32];
	__shared__ float tileB[32][32];

	int tx = threadIdx.x,ty = threadIdx.y;
	int n = blockIdx.x * blockDim.x + threadIdx.x;
	int m = blockIdx.y * blockDim.y + threadIdx.y;
	if(n >=N || m>=M)
	{
		return;
	}

	float psum = 0.0;
	for(int idx_tile = 0; idx_tile < K/32; idx_tile++)
	{
		tileA[ty][tx] = A[m*K + tx + idx_tile*32];
		tileB[ty][tx] = B[(idx_tile * 32 + ty) * K + n];
		__syncthreads();
		for(int k = 0; k < 32; k++)
		{
			psum += tileA[ty][k] * tileB[k][tx];
		}
		__syncthreads();
	}
	C[n + m*N] = psum;
}