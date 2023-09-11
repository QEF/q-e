#include <hip/hip_runtime.h>
#include <stdio.h>

__global__ void scalar_multiply_kernel(double *dev_ptr, double val, int s) 
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;
        int stride = gridDim.x * blockDim.x;
	for (; i < s; i += stride) {
		dev_ptr[i] *= val;
	}
}

extern "C" { 
  void c_scalar_multiply_(int s, double val, double *dev_ptr, hipStream_t stream)
  {
   
    const int blocks = 32;
    const int threads = 256;
    hipError_t err;

    //for (int i=0; i<s; ++i)
    //  printf("before --> %f\n",dev_ptr[i]);

    hipLaunchKernelGGL(scalar_multiply_kernel, dim3(blocks), dim3(threads), 0, stream, dev_ptr, val, s );
    //err = hipDeviceSynchronize();

    //for (int i=0; i<s; ++i)
    //  printf("after --> %f\n",dev_ptr[i]);

  }
}
