#include <stdio.h>
#include <stdlib.h>
#include "hostFE.h"
#include "helper.h"

void hostFE(int filterWidth, float *filter, int imageHeight, int imageWidth,
            float *inputImage, float *outputImage, cl_device_id *device,
            cl_context *context, cl_program *program)
{
    cl_int status;
    int filterSize = filterWidth * filterWidth * sizeof(float);
    int imgSize = imageHeight * imageWidth * sizeof(float);
    size_t localws[2] = {8, 8};
    size_t globalws[2] = {imageHeight, imageWidth};

    cl_command_queue queue;
    queue = clCreateCommandQueue(*context, *device, 0, NULL);

    cl_mem d_filterw = clCreateBuffer(*context, CL_MEM_USE_HOST_PTR, sizeof(int), &filterWidth, NULL);
    cl_mem d_filter = clCreateBuffer(*context, CL_MEM_USE_HOST_PTR, filterSize, filter, NULL);
    cl_mem d_in = clCreateBuffer(*context, CL_MEM_USE_HOST_PTR, imgSize, inputImage, NULL);
    cl_mem d_out = clCreateBuffer(*context, CL_MEM_WRITE_ONLY, imgSize, NULL, NULL);
    
    cl_kernel kernel = clCreateKernel(*program, "convolution", NULL);
    clSetKernelArg(kernel, 0, sizeof(cl_mem), (void*) &d_filterw);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), (void*) &d_filter);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), (void*) &d_in);
    clSetKernelArg(kernel, 3, sizeof(cl_mem), (void*) &d_out);

    clEnqueueNDRangeKernel(queue, kernel, 2, NULL, globalws, localws, 0, NULL, NULL);
    clEnqueueReadBuffer(queue, d_out, CL_TRUE, 0, imgSize, (void*)outputImage, 0, NULL, NULL);

}