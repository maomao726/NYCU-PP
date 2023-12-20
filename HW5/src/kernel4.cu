#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#define TILE_WIDTH 8

__global__ void mandelKernel( int* d_data,
                              int width,
                              float lowerX, float lowerY,
                              float stepX, float stepY,
                              int max_iteration, int pitch)
{
    // To avoid error caused by the floating number, use the following pseudo code
    //
    // float x = lowerX + thisX * stepX;
    // float y = lowerY + thisY * stepY;
    int thisX = blockIdx.x * blockDim.x + threadIdx.x;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;
    float c_re = lowerX + thisX * stepX;
    float c_im = lowerY + thisY * stepY;

    float z_re = c_re, z_im = c_im;
    int i;
    for(i = 0; i < max_iteration; i++)
    {
        if (z_re * z_re + z_im * z_im > 4.f)
            break;

        float new_re = z_re * z_re - z_im * z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }

    //write i back
    int* target = (int*)((char*)d_data + thisY * pitch) + thisX;
    *target = i;
}

// Host front-end function that allocates the memory and launches the GPU kernel
void hostFE (float upperX, float upperY, float lowerX, float lowerY, int* img, int resX, int resY, int maxIterations)
{
    float stepX = (upperX - lowerX) / resX;
    float stepY = (upperY - lowerY) / resY;

    int img_size = resX * resY * sizeof(int);
    int* h_data;
    int* d_data;
    size_t pitch;
    
    //kernel config/invoke
    cudaHostAlloc((void**) &h_data, img_size, cudaHostAllocDefault);
    cudaMallocPitch((void**) &d_data, &pitch, resX * sizeof(int), resY);
    dim3 dimGrid(resX / TILE_WIDTH, resY / TILE_WIDTH);
    dim3 dimBlock(TILE_WIDTH, TILE_WIDTH);
    mandelKernel<<<dimGrid, dimBlock>>>(d_data, resX, lowerX, lowerY, stepX, stepY, maxIterations, pitch);

    cudaMemcpy2D(h_data, resX * sizeof(int), d_data, pitch, resX * sizeof(int), resY, cudaMemcpyDeviceToHost);
    memcpy(img, h_data, img_size);

    cudaFree(d_data);
    cudaFreeHost(h_data);
}
